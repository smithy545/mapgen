//
// Created by Philip Smith on 4/17/2021.
//

#include <mapgen/terrain.h>

#include <BulletCollision/NarrowPhaseCollision/btRaycastCallback.h>
#include <cassert>
#include <engine/InstanceList.h>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/matrix_decompose.hpp>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>


using namespace engine;
using namespace utils::math;

namespace mapgen {
    Terrain::Terrain(unsigned int num_regions,
                     unsigned int width,
                     unsigned int height,
                     unsigned int num_tectonic_plates,
                     float wind_angle)
                     : Terrain(generate_sites(num_regions, width, height), num_tectonic_plates, wind_angle) {}

    Terrain::Terrain(std::vector<double> coords,
                     unsigned int num_tectonic_plates,
                     float wind_angle) : m_topology(coords) {
        assert(num_tectonic_plates > 0);


        for(const auto& region_topology: m_topology.get_regions())
            m_regions.emplace_back(region_topology);

        // assign ocean tiles
        assign_ocean_from_hull(2);

        // assign sea levels (sea level = path length to nearest ocean tiles)
        std::unordered_set<std::size_t> added;
        std::unordered_set<std::size_t> layer;
        for (auto index: m_oceans) {
            layer.insert(index);
            added.insert(index);
        }
        auto level = 0;
        do {
            std::unordered_set<std::size_t> next_layer;
            for (auto i: layer) {
                auto& region = m_regions.at(i);
                const auto& topology = region.get_topology();
                region.set_sea_level(level);
                for (auto n: topology.get_neighborhood()) {
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
        std::vector<std::size_t> tectonic_plates;
        while (tectonic_plates.size() < num_tectonic_plates) {
            auto index = rand() % m_regions.size();
            auto& region = m_regions[index];
            if(region.get_plate() == -1) {
                region.set_plate(tectonic_plates.size());
                tectonic_plates.push_back(index);
            }
        }

        // assign terrain regions to tectonic plates and make mountains along plate boundaries
        auto num_added = tectonic_plates.size();
        int mountain_size = 2;
        while(num_added < m_regions.size()) {
            std::vector<std::size_t> next;
            for (auto index : tectonic_plates) {
                auto& region = m_regions[index];
                const auto& topology = region.get_topology();
                for (auto n: topology.get_neighborhood()) {
                    auto& neighbor = m_regions[n];
                    if (neighbor.get_plate() == -1) {
                        neighbor.set_plate(region.get_plate());
                        num_added++;
                        next.push_back(n);
                    } else if(!m_oceans.contains(n) // tile can't be mountain and ocean
                    && neighbor.get_plate() != region.get_plate() // mountains only on tectonic plate boundaries
                    && neighbor.get_sea_level() > mountain_size + 2) // no mountains on beach
                        m_mountains.insert(n);
                }
            }
            tectonic_plates = next;
        }

        // assign elevations based on distance from mountain or distance from ocean (whichever is closer)
        for (auto i = 0; i < m_regions.size(); i++) {
            auto& region = m_regions[i];
            if (m_oceans.contains(i)) {
                region.set_elevation(0.0);
            } else if (m_mountains.contains(i)) {
                region.set_elevation(1.0);
            } else {
                int path_length;
                auto mountain_region = m_regions[find_nearest_mountain_face(i, path_length)];
                auto mountain_site = mountain_region.get_topology().get_site();
                auto dx = glm::cos(wind_angle);
                auto dy = glm::sin(wind_angle);
                auto site = region.get_topology().get_site();
                auto projection = dx*site.x + dy*site.y;
                auto mountain_projection = dx*mountain_site.x + dy*mountain_site.y;
                if(path_length > mountain_size) {
                    region.set_elevation((0.2f * region.get_sea_level()) / level);
                    region.set_moisture(0.8f);
                } else {
                    auto d = (10.0 * (path_length - 1)) / mountain_size;
                    region.set_elevation(0.8 / glm::log(d + glm::exp(1)));
                    if(projection > mountain_projection)
                        region.set_moisture(1.0f * (projection - mountain_projection)/(mountain_size*projection));
                    else
                        region.set_moisture(0.8f);
                }
            }
        }
    }

    entt::entity Terrain::register_mesh(entt::registry &registry) {
        Mesh mesh;
        std::unordered_map<std::string, std::size_t> added_points;
        for (const auto& region: m_regions) {
            auto index = region.get_index();
            const auto& topology = region.get_topology();
            auto ocean_color = glm::vec3(0, 0, 1);
            auto desert_color = glm::vec3(251, 225, 182)/255.f;
            auto grassland_color = glm::vec3(83, 91, 41)/255.f;
            auto mountain_color = glm::vec3(.6, .6, .6);
            auto mountain_top_color = glm::vec3(1, 1, 1);
            auto color = ocean_color;
            auto elevation = region.get_elevation();
            if (!m_oceans.contains(index)) {
                elevation *= MOUNTAIN_HEIGHT;
                if(region.get_sea_level() == 1 || region.get_moisture() < .3)
                    color = desert_color;
                else
                    color = grassland_color * region.get_moisture();
                if(elevation > .4f * MOUNTAIN_HEIGHT)
                    color = mountain_color * elevation*.5f;
                if (elevation > .6f * MOUNTAIN_HEIGHT)
                    color = mountain_color;
                if (elevation >= .9f * MOUNTAIN_HEIGHT)
                    color = mountain_top_color;
            }
            auto site = topology.get_site();
            mesh.vertices.emplace_back(site.x, elevation, site.y);
            mesh.colors.push_back(color);

            std::map<double, std::size_t> neighbors;
            for (auto n: topology.get_neighborhood()) {
                const auto& topology = m_topology.get_region_topology(n);
                auto neighbor_site = topology.get_site();
                neighbors[std::atan2(
                        neighbor_site.y - site.y,
                        neighbor_site.x - site.x)] = n;
            }
            auto last_angle = glm::two_pi<double>();
            for (auto[angle, n1]: neighbors) {
                if (last_angle != glm::two_pi<double>()) {
                    auto n2 = neighbors[last_angle];
                    mesh.indices.push_back(region.get_index());
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

        return m_entity;
    }

    entt::entity Terrain::register_wireframe_mesh(entt::registry &registry) {
        Mesh mesh;
        std::unordered_map<std::string, std::size_t> added_points;
        for (auto& region: m_regions) {
            const auto& topology = region.get_topology();
            auto elevation = region.get_elevation();
            auto site = topology.get_site();
            if (!m_oceans.contains(region.get_index()))
                elevation *= MOUNTAIN_HEIGHT;
            mesh.vertices.emplace_back(site.x, elevation, site.y);
            if(m_rivers.contains(mesh.vertices.size() - 1))
                mesh.colors.emplace_back(.1, .2, 1.0);
            else
                mesh.colors.emplace_back(0,0,0);
            for (auto n: topology.get_neighborhood()) {
                if(n > region.get_index()) {
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

    void Terrain::assign_ocean_from_hull(int depth) {
        for(auto index: m_topology.get_hull_regions()) {
            apply_callback_breadth_first(
                    [&](RegionTerrain& region, int depth) {
                        m_oceans.insert(region.get_index());
                        region.set_moisture(1.0);
                    },
                    m_regions[index],
                    depth,
                    [&](RegionTerrain& region) {
                        return !m_oceans.contains(region.get_index());
                    });
        }
    }

    void Terrain::apply_callback_breadth_first(const std::function<void(RegionTerrain&, std::size_t)>& callback,
                                               RegionTerrain& start,
                                               int max_depth,
                                               const std::function<bool(RegionTerrain&)>& predicate,
                                               std::unordered_set<std::size_t> visited) {
        int depth = 0;
        std::unordered_set<std::size_t> current{start.get_index()};
        while(depth < max_depth) {
            std::unordered_set<std::size_t> next;
            for(auto i: current) {
                auto& region = m_regions[i];
                const auto& topology = region.get_topology();
                callback(region, i);
                visited.insert(i);
                for (auto neighbor: topology.get_neighborhood()) {
                    if (!visited.contains(neighbor) && !next.contains(neighbor) && predicate(m_regions[neighbor]))
                        next.insert(neighbor);
                }
            }
            current = next;
            depth++;
        }
    }

    void Terrain::apply_callback_sweep_line(const std::function<void(RegionTerrain&, float)>& callback,
                                            float direction_in_radians,
                                            const std::function<bool(RegionTerrain&, float)>& predicate) {
        auto dx = glm::cos(direction_in_radians);
        auto dy = glm::sin(direction_in_radians);
        std::map<float, RegionTerrain&> ordered_regions;
        for(auto &region: m_regions) {
            const auto& topology = region.get_topology();
            auto site = topology.get_site();
            auto projection = dx*site.x + dy*site.y;
            ordered_regions.insert({projection, region});
        }
        for(auto[distance, region]: ordered_regions) {
            if(predicate(region, distance))
                callback(region, distance);
        }
    }

    std::size_t Terrain::find_nearest_mountain_face(std::size_t index, int& path_length) {
        assert(!m_mountains.empty());
        if(m_mountains.contains(index))
            return index;
        std::unordered_set<std::size_t> layer;
        std::unordered_set<std::size_t> searched{index};
        auto min_distance = std::numeric_limits<double>::infinity();
        auto base = m_topology.get_region_topology(index);
        for(auto n: base.get_neighborhood())
            layer.insert(n);
        auto found{index};
        int level{0};
        while(found == index && searched.size() < m_regions.size()) {
            level++;
            std::unordered_set<std::size_t> next;
            for (auto i: layer) {
                const auto& topology = m_topology.get_region_topology(i);
                if (m_mountains.contains(i)) {
                    auto distance = glm::distance(topology.get_site(), base.get_site());
                    if (distance < min_distance) {
                        found = i;
                        min_distance = distance;
                    }
                } else {
                    for(auto n: topology.get_neighborhood()) {
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

    std::vector<double> Terrain::generate_sites(unsigned int num_regions, unsigned int width, unsigned int height) {
        assert(num_regions > 0);
        assert(width > 0);
        assert(height > 0);
        std::random_device rd;
        std::vector<double> coords;
        std::unordered_set<glm::vec2> sites;
        while(sites.size() < num_regions) {
            glm::vec2 v;
            do {
                v = glm::vec2(rd() % width, rd() % height);
            } while (sites.contains(v));
            sites.insert(v);
            coords.push_back(v.x);
            coords.push_back(v.y);
        }
        return coords;
    }

    RegionTerrain::RegionTerrain(const RegionTopology& neighbor_topology) : m_topology(neighbor_topology) {}

    std::size_t RegionTerrain::get_index() const {
        return m_topology.get_index();
    }
} // namespace mapgen
