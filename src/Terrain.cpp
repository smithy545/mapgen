//
// Created by Philip Smith on 4/17/2021.
//

#include <mapgen/Terrain.h>

#include <BulletCollision/NarrowPhaseCollision/btRaycastCallback.h>
#include <cassert>
#include <engine/InstanceList.h>
#include <engine/Mesh.h>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/matrix_decompose.hpp>
#include <limits>
#include <map>
#include <mapgen/DelaunatorAlgorithm.h>
#include <random>
#include <unordered_map>
#include <utils/math_util.h>
#include <utility>
#include <vector>


using namespace engine;
using namespace utils::math;

namespace mapgen {
    Terrain::Terrain(unsigned int num_sites, int width, int height, int num_tectonic_plates, bool centered) {
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
                if (centered)
                    v = glm::vec2(rd() % width - width / 2, rd() % height - height / 2);
                else
                    v = glm::vec2(rd() % width, rd() % height);
            } while (site_hashes.contains(Diagram::site_key(v)));
            site_hashes.insert(Diagram::site_key(v));
            coords.push_back(v.x);
            coords.push_back(v.y);
        }
        m_base = DelaunatorAlgorithm::construct_clamped_voroni_diagram(coords, 0, 0, width, height);
        // apply lloyd relaxation for prettiness
        m_base = m_base.relax();
        m_dual = m_base.dual();

        // randomly choose tectonic plate origins from faces
        while(tectonic_plates.size() < num_tectonic_plates)
            tectonic_plates.insert(rand() % m_base.get_faces().size());

        // assign ocean tiles
        for (auto index: m_base.get_hull())
            assign_ocean(index, 2);
    }

    void Terrain::register_voroni_debug_mesh(entt::registry &registry) {
        Mesh edge_mesh, face_mesh, site_mesh;
        for(auto edge: m_base.get_edges()) {
            auto i = edge_mesh.vertices.size();
            edge_mesh.vertices.emplace_back(edge.first.x, -1, edge.first.y);
            edge_mesh.vertices.emplace_back(edge.second.x, -1, edge.second.y);
            edge_mesh.colors.emplace_back(0,0,1);
            edge_mesh.colors.emplace_back(0,0,1);
            edge_mesh.indices.push_back(i++);
            edge_mesh.indices.push_back(i);
        }
        auto e1 = registry.create();
        registry.emplace_or_replace<Mesh>(e1, edge_mesh);
        registry.patch<InstanceList>(e1, [](auto &instance_list) {
            instance_list.set_instances(std::vector<glm::mat4>{glm::mat4(1)});
            instance_list.render_strategy = GL_LINES;
        });
        for(const auto &face: m_base.get_faces()) {
            for(auto [n, e]: face.neighboring_edges) {
                auto edge = m_base.get_edges()[e];
                auto i = face_mesh.vertices.size();
                face_mesh.vertices.emplace_back(edge.first.x, 0, edge.first.y);
                face_mesh.vertices.emplace_back(edge.second.x, 0, edge.second.y);
                face_mesh.colors.emplace_back(1, 1, 0);
                face_mesh.colors.emplace_back(1, 1, 0);
                face_mesh.indices.push_back(i++);
                face_mesh.indices.push_back(i);
            }
            auto i = site_mesh.vertices.size();
            site_mesh.vertices.emplace_back(face.site.x, 0, face.site.y);
            site_mesh.colors.emplace_back(1,0,0);
            site_mesh.indices.push_back(i);
        }
        auto e2 = registry.create();
        registry.emplace_or_replace<Mesh>(e2, face_mesh);
        registry.patch<InstanceList>(e2, [](auto &instance_list) {
            instance_list.set_instances(std::vector<glm::mat4>{glm::mat4(1)});
            instance_list.render_strategy = GL_LINES;
        });
        auto e3 = registry.create();
        registry.emplace_or_replace<Mesh>(e3, site_mesh);
        registry.patch<InstanceList>(e3, [](auto &instance_list) {
            instance_list.set_instances(std::vector<glm::mat4>{glm::mat4(1)});
            instance_list.render_strategy = GL_POINTS;
        });
    }

    void Terrain::register_terrain_mesh(entt::registry &registry) {
        // generate terrain regions
        for (unsigned int i = 0; i < m_base.get_faces().size(); i++) {
            regions.insert({i, TerrainRegion{i}});
        }

        // assign layers where each layer's invariant is the distance to the closest ocean site
        std::unordered_set<unsigned int> added;
        std::unordered_set<unsigned int> layer;
        for (auto index: ocean) {
            for (auto[n, e]: m_base.get_faces()[index].neighboring_edges) {
                layer.insert(n);
                added.insert(n);
            }
        }

        auto level = 0.0;
        do {
            std::unordered_set<unsigned int> next_layer;
            for (auto i: layer) {
                regions[i].elevation = level;
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

        // create tectonic plates
        added.clear();
        std::vector<unsigned int> ordered;
        ordered.insert(ordered.end(), tectonic_plates.begin(), tectonic_plates.end());
        std::map<unsigned int, std::unordered_set<unsigned int>> faces;
        for(auto plate: ordered) {
            faces.insert({plate, {plate}});
            added.insert(plate);
        }

        // generate and assign terrain regions to tectonic plates and make mountains along plate boundaries
        while(regions.size() != m_base.get_faces().size()) {
            std::map<unsigned int, std::unordered_set<unsigned int>> next;
            for(auto plate : ordered) {
                next[plate] = {};
                for(auto face_index: faces[plate]) {
                    auto &face = m_base.get_faces()[face_index];
                    regions[face_index].plate = plate;
                    for (auto[n, e]: face.neighboring_edges) {
                        if (added.contains(n)) {
                            if(regions[n].plate != plate)
                                mountains.insert(n);
                        } else {
                            bool found{false};
                            for(auto [k, v]: next) {
                                if(k != plate && v.contains(n)) {
                                    found = true;
                                    break;
                                }
                            }
                            if(!found)
                                next[plate].insert(n);
                            added.insert(n);
                        }
                    }
                }
            }
            faces = next;
        }

        for (auto [index, region]: regions) {
            auto& face = m_base.get_faces()[index];
            if (!ocean.contains(index) && !mountains.contains(index)) {
                auto nearest_index = find_nearest_mountain_face_recursive(index);
                if (nearest_index != index) {
                    auto nearest = m_base.get_faces()[nearest_index];
                    auto d = glm::distance(nearest.site, face.site);
                    // scale height inverse logarithmically to distance from mountains (range 0-1 when shifting x by e)
                    if (d > 5.0)
                        region.elevation *= 2.0 / glm::log(d + glm::exp(1));
                    else
                        region.elevation = regions[nearest_index].elevation; // turn into plateau if close enough
                }
            }
        }

        // generate meshes
        Mesh mesh;
        std::unordered_map<std::string, unsigned int> added_points;
        for (auto [index, region]: regions) {
            if (!ocean.contains(index))
                region.elevation *= 200;
            auto ocean_color = glm::vec3(0, 0, 1);
            auto land_color = glm::vec3(.5, .3, .1);
            auto mountain_color = glm::vec3(.6, .6, .6);
            auto mountain_top_color = glm::vec3(1, 1, 1);
            auto color = ocean_color;
            if (region.elevation > 0)
                color = land_color + glm::vec3(.2, .4, .4) / (float)region.elevation;
            if (region.elevation > 80)
                color = mountain_color;
            if (region.elevation > 200)
                color = mountain_top_color;
            auto& face = m_base.get_faces()[index];
            mesh.vertices.emplace_back(face.site.x, region.elevation, face.site.y);
            mesh.colors.push_back(color);

            std::map<double, unsigned int> neighbors;
            for (auto[n, e]: face.neighboring_edges) {
                auto& neighbor = m_base.get_faces()[n];
                neighbors[std::atan2(
                        neighbor.site.y - face.site.y,
                        neighbor.site.x - face.site.x)] = n;
            }
            auto last_angle = glm::two_pi<double>();
            for (auto[angle, n1]: neighbors) {
                if (last_angle != glm::two_pi<double>()) {
                    auto n2 = neighbors[last_angle];
                    mesh.indices.push_back(index);
                    mesh.indices.push_back(n1);
                    mesh.indices.push_back(n2);
                }
                last_angle = angle;
            }
        }

        m_entity = registry.create();
        // register render data
        registry.emplace_or_replace<Mesh>(m_entity, mesh);
        registry.patch<InstanceList>(m_entity, [](auto &instance_list) {
            instance_list.set_instances(std::vector<glm::mat4>{glm::mat4(1)});
        });
        // register physics data
        m_mesh = std::make_shared<btTriangleMesh>();
        for(auto i = 0; i < mesh.indices.size(); i += 3) {
            auto v1 = glm2bt(mesh.vertices[mesh.indices[i]]);
            auto v2 = glm2bt(mesh.vertices[mesh.indices[i + 1]]);
            auto v3 = glm2bt(mesh.vertices[mesh.indices[i + 2]]);
            m_mesh->addTriangle(v1, v2, v3);
        }
        m_shape = std::make_shared<btBvhTriangleMeshShape>(m_mesh.get(), true);
        m_body = std::make_shared<btCollisionObject>();
        m_body->setCollisionShape(m_shape.get());
    }

    void Terrain::assign_ocean(unsigned int start, int neighbor_depth) {
        if (neighbor_depth == 0)
            return;
        ocean.insert(start);
        for (auto[neighbor, edge]: m_base.get_faces()[start].neighboring_edges)
            assign_ocean(neighbor, neighbor_depth - 1);
    }

    unsigned int Terrain::find_nearest_mountain_face(unsigned int index) {
        auto distance = std::numeric_limits<double>::infinity();
        auto face = m_base.get_faces()[index];
        auto found = index;
        for (auto i: mountains) {
            if (index != i) {
                auto mountain_face = m_base.get_faces()[i];
                auto d = glm::distance(face.site, mountain_face.site);
                if (found == index || d < distance) {
                    found = i;
                    distance = d;
                }
            }
        }
        return found;
    }

    unsigned int Terrain::find_nearest_mountain_face_recursive(unsigned int index) {
        if(mountains.empty())
            return index;
        auto face = m_base.get_faces()[index];
        auto neighbors = face.neighboring_edges;
        std::unordered_set<unsigned int> next{};
        std::unordered_set<unsigned int> searched{index};
        auto distance = std::numeric_limits<double>::infinity();
        auto found = index;
        for (auto[n, e]: neighbors) {
            if (mountains.contains(n) && glm::distance(face.site, m_base.get_faces()[n].site) < distance) {
                found = n;
                distance = glm::distance(face.site, m_base.get_faces()[n].site);
            }
            searched.insert(n);
            for (auto[n2, e2]: m_base.get_faces()[n].neighboring_edges) {
                if (!searched.contains(n2))
                    next.insert(n2);
            }
        }
        if (found == index)
            return find_mountain_kernel(index, next, searched);
        return found;
    }

    unsigned int Terrain::find_mountain_kernel(unsigned int index, const std::unordered_set<unsigned int> &to_search,
                                               std::unordered_set<unsigned int> searched) {
        if (searched.size() == m_base.get_faces().size())
            return index;
        std::unordered_set<unsigned int> next{};
        auto distance = std::numeric_limits<double>::infinity();
        auto found = index;
        auto base = m_base.get_faces()[index];
        for (auto i: to_search) {
            auto face = m_base.get_faces()[i];
            if (mountains.contains(i) && glm::distance(face.site, base.site) < distance) {
                found = i;
                distance = glm::distance(face.site, base.site);
            }
            searched.insert(i);
            for (auto[n, e]: face.neighboring_edges) {
                if (!searched.contains(n))
                    next.insert(n);
            }
        }
        if (found == index)
            find_mountain_kernel(index, next, searched);
        return found;
    }

    std::vector<glm::vec3> Terrain::get_mouse_terrain_collision(float x, float y, const RenderContext& context) const {
        if(m_mesh == nullptr || m_shape == nullptr) {
            return {};
        }
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
