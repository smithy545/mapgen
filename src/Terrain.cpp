//
// Created by Philip Smith on 4/17/2021.
//

#include <mapgen/Terrain.h>

#include <BulletCollision/NarrowPhaseCollision/btRaycastCallback.h>
#include <cassert>
#include <engine/InstanceList.h>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/matrix_decompose.hpp>
#include <glm/gtx/projection.hpp>
#include <limits>
#include <map>
#include <mapgen/DelaunatorAlgorithm.h>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utils/math_util.h>
#include <utility>
#include <vector>


using namespace engine;
using namespace utils::math;

namespace mapgen {
    Terrain::Terrain(unsigned int num_sites, int width, int height, int num_tectonic_plates, glm::vec2 wind_direction) {
        assert(num_sites > 0);
        assert(width > 0);
        assert(height > 0);
        assert(num_tectonic_plates > 0);
        std::random_device rd;
        std::vector<double> coords;
        std::unordered_set<std::string> site_hashes;
        for (int i = 0; i < num_sites; i++) {
            glm::vec2 v;
            do {
                v = glm::vec2(rd() % width, rd() % height);
            } while (site_hashes.contains(Diagram::site_key(v)));
            site_hashes.insert(Diagram::site_key(v));
            coords.push_back(v.x);
            coords.push_back(v.y);
        }
        m_base = DelaunatorAlgorithm::construct_voroni_diagram(coords);
        m_dual = DelaunatorAlgorithm::construct_delauney_diagram(coords);

        // generate terrain regions over geometry
        for (const auto& face: m_base.get_faces())
            m_regions.push_back(TerrainRegion{face});

        // assign ocean tiles
        assign_ocean_from_hull(2);

        // assign sea levels (sea level = distance to nearest ocean tiles)
        std::unordered_set<unsigned int> added;
        std::unordered_set<unsigned int> layer;
        for (auto index: m_oceans) {
            layer.insert(index);
            added.insert(index);
        }
        auto level = 0;
        do {
            std::unordered_set<unsigned int> next_layer;
            for (auto i: layer) {
                m_regions[i].sea_level = level;
                for (auto[n, e]: m_base.get_faces()[i].neighboring_edges) {
                    if (!added.contains(n)) {
                        next_layer.insert(n);
                        added.insert(n);
                    }
                }
            }
            layer = next_layer;
            level += 1;
        } while (!layer.empty());

        // randomly choose tectonic plate origins from faces
        std::vector<unsigned int> tectonic_plates;
        while (tectonic_plates.size() < num_tectonic_plates) {
            auto index = rand() % m_base.get_faces().size();
            if(m_regions[index].plate_id == -1) {
                m_regions[index].plate_id = tectonic_plates.size();
                tectonic_plates.push_back(index);
            }
        }

        // assign terrain regions to tectonic plates and make mountains along plate boundaries
        int num_added = tectonic_plates.size();
        int mountain_size = 2;
        while(num_added < m_regions.size()) {
            std::vector<unsigned int> next;
            for (auto index : tectonic_plates) {
                auto& region = m_regions[index];
                for (auto[n, e]: region.face.neighboring_edges) {
                    if (m_regions[n].plate_id == -1) {
                        m_regions[n].plate_id = region.plate_id;
                        num_added++;
                        next.push_back(n);
                    } else if(!m_oceans.contains(n) // no mountain/ocean combo tiles
                    && m_regions[n].plate_id != region.plate_id // mountains only on tectonic plate boundaries
                    && m_regions[n].sea_level > mountain_size + 2) // no mountains on beach
                        m_mountains.insert(n);
                }
            }
            tectonic_plates = next;
        }

        // assign elevations based on distance from mountain or distance from ocean (whichever is closer)
        for (auto i = 0; i < m_regions.size(); i++) {
            auto &region = m_regions[i];
            if (m_oceans.contains(i)) {
                region.elevation = 0.0;
            } else if (m_mountains.contains(i)) {
                region.elevation = 1.0;
            } else {
                int path_length;
                find_nearest_mountain_face(i, path_length);
                if(path_length > mountain_size)
                    region.elevation = (0.2 * region.sea_level) / level;
                else {
                    auto d = (10.0 * (path_length - 1)) / mountain_size;
                    region.elevation = 0.8 / glm::log(d + glm::exp(1));
                }
                region.moisture = 1.0 - region.elevation;
            }
        }

        // rivers
        std::unordered_map<unsigned int, unsigned int> nearest_mountain;
        for(auto index: m_mountains) {
            auto& region = m_regions[index];
            auto wind_projection = glm::length(glm::proj(region.face.site, wind_direction));
            apply_callback_breadth_first(
                    [&](unsigned int index) {

                    },
                    index,
                    2,
                    [&](unsigned int index) {
                        return !m_mountains.contains(index) && !m_oceans.contains(index);
                    });
        }
    }

    void Terrain::register_mesh(entt::registry &registry) {
        Mesh mesh;
        std::unordered_map<std::string, unsigned int> added_points;
        for (auto& region: m_regions) {
            auto ocean_color = glm::vec3(0, 0, 1);
            auto desert_color = glm::vec3(251./255., 225./255., 182./255.);
            auto grassland_color = glm::vec3(83./255., 91./255., 41./255.);
            auto mountain_color = glm::vec3(.6, .6, .6);
            auto mountain_top_color = glm::vec3(1, 1, 1);
            auto color = ocean_color;
            auto elevation = region.elevation;
            if (!m_oceans.contains(region.face.id)) {
                elevation *= 200;
                if(region.sea_level == 1 || region.moisture < .3)
                    color = desert_color;
                else
                    color = grassland_color * region.moisture;
                if(elevation > 80)
                    color = mountain_color * elevation*.5f;
                if (elevation > 120)
                    color = mountain_color;
                if (elevation >= 180)
                    color = mountain_top_color;
            }
            mesh.vertices.emplace_back(region.face.site.x, elevation, region.face.site.y);
            mesh.colors.push_back(color);

            std::map<double, unsigned int> neighbors;
            for (auto[n, e]: region.face.neighboring_edges) {
                auto& neighbor = m_base.get_faces()[n];
                neighbors[std::atan2(
                        neighbor.site.y - region.face.site.y,
                        neighbor.site.x - region.face.site.x)] = n;
            }
            auto last_angle = glm::two_pi<double>();
            for (auto[angle, n1]: neighbors) {
                if (last_angle != glm::two_pi<double>()) {
                    auto n2 = neighbors[last_angle];
                    mesh.indices.push_back(region.face.id);
                    mesh.indices.push_back(n1);
                    mesh.indices.push_back(n2);
                }
                last_angle = angle;
            }
        }

        m_entity = registry.create();
        registry.emplace_or_replace<Mesh>(m_entity, mesh);
        registry.patch<InstanceList>(m_entity, [](auto &instance_list) {
            instance_list.set_instances(std::vector<glm::mat4>{glm::mat4(1)});
        });

        // register physics data
        m_mesh = std::make_shared<btTriangleMesh>();
        for (auto i = 0; i < mesh.indices.size(); i += 3) {
            auto v1 = glm2bt(mesh.vertices[mesh.indices[i]]);
            auto v2 = glm2bt(mesh.vertices[mesh.indices[i + 1]]);
            auto v3 = glm2bt(mesh.vertices[mesh.indices[i + 2]]);
            m_mesh->addTriangle(v1, v2, v3);
        }
        m_shape = std::make_shared<btBvhTriangleMeshShape>(m_mesh.get(), true);
        m_body = std::make_shared<btCollisionObject>();
        m_body->setCollisionShape(m_shape.get());
    }

    entt::entity Terrain::register_wireframe_mesh(entt::registry &registry) {
        Mesh mesh;
        std::unordered_map<std::string, unsigned int> added_points;
        for (auto& region: m_regions) {
            auto elevation = region.elevation;
            if (!m_oceans.contains(region.face.id))
                elevation *= 200;
            mesh.vertices.emplace_back(region.face.site.x, elevation, region.face.site.y);
            mesh.colors.emplace_back(0,0,0);
            for (auto[n, e]: region.face.neighboring_edges) {
                if(n > region.face.id) {
                    mesh.indices.push_back(mesh.vertices.size() - 1);
                    mesh.indices.push_back(n);
                }
            }
        }
        auto entity = registry.create();
        registry.emplace_or_replace<Mesh>(entity, mesh);
        registry.patch<InstanceList>(entity, [](auto &instance_list) {
            instance_list.set_instances(std::vector<glm::mat4>{glm::mat4(1)});
            instance_list.render_strategy = GL_LINES;
        });
        return entity;
    }

    entt::entity Terrain::register_voroni_mesh(entt::registry &registry) {
        Mesh face_mesh;
        for(const auto &face: m_base.get_faces()) {
            for(auto [n, e]: face.neighboring_edges) {
                if(n > face.id) {
                    auto edge = m_base.get_edges()[e];
                    auto i = face_mesh.vertices.size();
                    face_mesh.vertices.emplace_back(edge.first.x, 0, edge.first.y);
                    face_mesh.vertices.emplace_back(edge.second.x, 0, edge.second.y);
                    face_mesh.colors.emplace_back(1, 1, 0);
                    face_mesh.colors.emplace_back(1, 1, 0);
                    face_mesh.indices.push_back(i++);
                    face_mesh.indices.push_back(i);
                }
            }
        }
        auto entity = registry.create();
        registry.emplace_or_replace<Mesh>(entity, face_mesh);
        registry.patch<InstanceList>(entity, [](auto &instance_list) {
            instance_list.set_instances(std::vector<glm::mat4>{glm::mat4(1)});
            instance_list.render_strategy = GL_LINES;
        });
        return entity;
    }

    Mesh Terrain::get_region_mesh(unsigned int index) {
        assert(index < m_regions.size());
        Mesh mesh;
        glm::vec3 color{1,1,0};
        auto& region = m_regions[index];
        mesh.vertices.emplace_back(region.face.site.x, region.elevation, region.face.site.y);
        mesh.colors.push_back(color);

        std::map<double, unsigned int> neighbors;
        for (auto[n, e]: region.face.neighboring_edges) {
            auto& neighbor = m_base.get_faces()[n];
            neighbors[std::atan2(
                    neighbor.site.y - region.face.site.y,
                    neighbor.site.x - region.face.site.x)] = n;
        }
        bool first{true};
        for (auto[angle, n1]: neighbors) {
            auto& neighbor = m_base.get_faces()[n1];
            int i = mesh.vertices.size();
            mesh.vertices.emplace_back(neighbor.site.x, region.elevation, neighbor.site.y);
            mesh.colors.push_back(color);
            if (!first) {
                mesh.indices.push_back(0);
                mesh.indices.push_back(i);
                mesh.indices.push_back(i - 1);
            }
            first = false;
        }
        mesh.indices.push_back(0);
        mesh.indices.push_back(mesh.vertices.size() - 1);
        mesh.indices.push_back(1);
        return mesh;
    }

    void Terrain::assign_ocean_from_hull(int depth) {
        for(auto index: m_base.get_hull())
            apply_callback_breadth_first(
            [&](auto && index) {
                m_oceans.insert(index);
                m_regions[index].moisture = 1.0;
            },
            index,
            depth,
            [&](unsigned int index) {
                return !m_oceans.contains(index);
            });
    }

    void Terrain::apply_callback_breadth_first(const std::function<void(unsigned int)>& callback,
                                               unsigned int start,
                                               int max_depth,
                                               const std::function<bool(unsigned int)>& predicate,
                                               std::unordered_set<unsigned int> visited) {
        int depth = 0;
        std::unordered_set<unsigned int> next, current{start};
        while(depth < max_depth) {
            for(auto i: current) {
                callback(i);
                visited.insert(i);
                for (auto[neighbor, edge]: m_base.get_faces()[i].neighboring_edges) {
                    if (!visited.contains(neighbor) && !next.contains(neighbor) && predicate(neighbor))
                        next.insert(neighbor);
                }
            }
            current = next;
            depth++;
        }
    }

    void Terrain::apply_callback_breadth_first_recursively(const std::function<void(unsigned int)>& callback,
                                                           unsigned int start,
                                                           int max_depth,
                                                           const std::function<bool(unsigned int)>& predicate,
                                                           std::unordered_set<unsigned int> visited) {
        if (max_depth == 0)
            return;
        callback(start);
        visited.insert(start);
        for (auto[neighbor, edge]: m_base.get_faces()[start].neighboring_edges) {
            if(!visited.contains(neighbor) && predicate(neighbor))
                apply_callback_breadth_first_recursively(callback,
                                                         neighbor,
                                                         max_depth - 1,
                                                         predicate,
                                                         visited);
        }
    }

    unsigned int Terrain::find_nearest_mountain_face(unsigned int index, int& path_length) {
        assert(!m_mountains.empty());
        if(m_mountains.contains(index))
            return index;
        std::unordered_set<unsigned int> layer;
        std::unordered_set<unsigned int> searched{index};
        auto min_distance = std::numeric_limits<double>::infinity();
        auto base = m_base.get_faces()[index];
        for(auto [n, e]: base.neighboring_edges)
            layer.insert(n);
        unsigned int found{index};
        int level{0};
        while(found == index && searched.size() < m_regions.size()) {
            level++;
            std::unordered_set<unsigned int> next;
            for (auto i: layer) {
                auto face = m_base.get_faces()[i];
                if (m_mountains.contains(i)) {
                    auto distance = glm::distance(face.site, base.site);
                    if (distance < min_distance) {
                        found = i;
                        min_distance = distance;
                    }
                } else {
                    for(auto [n, e]: face.neighboring_edges) {
                        if(!searched.contains(n))
                            next.insert(n);
                    }
                }
                searched.insert(i);
            }
            layer = next;
        }
        path_length = level;
        return found;
    }

    glm::vec3 Terrain::get_mouse_terrain_collision(float x, float y, const RenderContext& context) const {
        if(m_mesh == nullptr || m_shape == nullptr)
            return {};
        auto cam = context.camera;
        auto position = cam->get_position();
        auto aspect = context.screen_width/context.screen_height;
        auto perspective = glm::perspective(context.fovy, aspect, context.z_near, context.z_far);
        glm::vec3 from{
                2.0 * (x - context.screen_width / 2) / context.screen_width,
                2.0 * (context.screen_height / 2 - y) / context.screen_height, 0};
        auto transformation = glm::inverse(perspective * cam->get_view()) * glm::translate(glm::mat4(1), from);
        glm::vec3 scale, skew;
        glm::vec4 persp;
        glm::quat rotation;
        glm::decompose(transformation, scale, rotation, from, skew, persp);
        auto to = from + (from - position) * context.z_far;

        btVector3 btFrom{from.x, from.y, from.z};
        btVector3 btTo{to.x, to.y, to.z};

        struct MyRaycastCallback : public btTriangleRaycastCallback {
        public:
            int m_index{-1};

            MyRaycastCallback(const btVector3& from, const btVector3& to)
                    : btTriangleRaycastCallback(from, to) {}

            btScalar reportHit(const btVector3& hitNormalLocal, btScalar hitFraction, int partId, int triangleIndex)
            override {
                m_index = triangleIndex;
                if(hitFraction < m_hitFraction)
                    return hitFraction;
                return m_hitFraction;
            }
        } res{btFrom, btTo};

        m_shape->performRaycast(&res, btFrom, btTo);
        if(res.m_index != -1)
            return from + res.m_hitFraction*(from - to);
        throw std::runtime_error("Mouse doesn't collide with terrain");
    }

    std::vector<glm::vec3> Terrain::get_mouse_terrain_colliding_triangle(float x,
                                                                         float y,
                                                                         const RenderContext& context) const {
        if(m_mesh == nullptr || m_shape == nullptr)
            return {};
        auto cam = context.camera;
        auto position = cam->get_position();
        auto aspect = context.screen_width/context.screen_height;
        auto perspective = glm::perspective(context.fovy, aspect, context.z_near, context.z_far);
        glm::vec3 from{
                2.0 * (x - context.screen_width / 2) / context.screen_width,
                2.0 * (context.screen_height / 2 - y) / context.screen_height, 0};
        auto transformation = glm::inverse(perspective * cam->get_view()) * glm::translate(glm::mat4(1), from);
        glm::vec3 scale, skew;
        glm::vec4 persp;
        glm::quat rotation;
        glm::decompose(transformation, scale, rotation, from, skew, persp);
        auto to = from + glm::normalize(from - position) * context.z_far;

        btVector3 bt_from{from.x, from.y, from.z};
        btVector3 bt_to{to.x, to.y, to.z};

        struct MyRaycastCallback : public btTriangleRaycastCallback {
        public:
            int m_index{-1};

            MyRaycastCallback(const btVector3& from, const btVector3& to)
                    : btTriangleRaycastCallback(from, to) {}

            btScalar reportHit(const btVector3& normal, btScalar fraction, int part, int index)
            override {
                m_index = index;
                if(fraction < m_hitFraction)
                    return fraction;
                return m_hitFraction;
            }
        } res{bt_from, bt_to};

        m_shape->performRaycast(&res, bt_from, bt_to);
        if(res.m_index != -1) {
            std::vector<glm::vec3> verts{};
            auto& tmesh = m_mesh->getIndexedMeshArray()[0];
            auto* gfxbase = (unsigned int*)(tmesh.m_triangleIndexBase+res.m_index*tmesh.m_triangleIndexStride);
            for (int j = 0; j < 3; j++) {
                int graphicsindex = tmesh.m_indexType==PHY_SHORT?((unsigned short*)gfxbase)[j]:gfxbase[j];
                auto* graphicsbase = (float*)(tmesh.m_vertexBase+graphicsindex*tmesh.m_vertexStride);
                btVector3 v(graphicsbase[0]*m_mesh->getScaling().getX(),
                            graphicsbase[1]*m_mesh->getScaling().getY(),
                            graphicsbase[2]*m_mesh->getScaling().getZ());
                verts.emplace_back(v.x(), v.y(), v.z());
            }

            return verts;
        }
        return {};
    }

} // namespace mapgen
