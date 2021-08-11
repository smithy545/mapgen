//
// Created by Philip Smith on 4/17/2021.
//

#include <mapgen/Terrain.h>

#include <BulletCollision/NarrowPhaseCollision/btRaycastCallback.h>
#include <cassert>
#include <engine/InstanceList.h>
#include <glm/glm.hpp>
#include <glm/gtx/matrix_decompose.hpp>
#include <limits>
#include <random>
#include <stdexcept>
#include <utils/math_util.h>
#include <vector>


using namespace engine;
using namespace utils::math;

namespace mapgen {
    Terrain::Terrain(unsigned int num_regions,
                     unsigned int width,
                     unsigned int height,
                     unsigned int num_tectonic_plates,
                     float wind_direction)
                     : Terrain(generate_sites(num_regions, width, height), num_tectonic_plates, wind_direction) {}

    Terrain::Terrain(std::vector<double> site_coords, unsigned int num_tectonic_plates, float wind_direction)
    : m_delauney(site_coords) {
        assert(num_tectonic_plates > 0);

        m_voroni_sites = std::move(site_coords);

        // TODO: Lock region index to half of coord index? Wastes memory when regions become deletable or mutable
        for(std::size_t i = 0; i < m_voroni_sites.size(); i += 2)
            m_regions.push_back(Region{i});

        std::unordered_map<glm::vec2, std::unordered_set<std::size_t>> edge_to_sites;
        std::unordered_map<std::size_t, glm::vec2> site_to_edge;
        for (std::size_t i = 0; i < m_delauney.triangles.size(); i += 3) {
            auto i0 = 2 * m_delauney.triangles[i];
            auto i1 = 2 * m_delauney.triangles[i + 1];
            auto i2 = 2 * m_delauney.triangles[i + 2];
            glm::vec2 v0(m_voroni_sites[i0], m_voroni_sites[i0 + 1]);
            glm::vec2 v1 (m_voroni_sites[i1], m_voroni_sites[i1 + 1]);
            glm::vec2 v2(m_voroni_sites[i2], m_voroni_sites[i2 + 1]);
            auto cc = utils::math::compute_triangle_circumcenter(v0, v1, v2);
            site_to_edge[i] = cc;
            site_to_edge[i + 1] = cc;
            site_to_edge[i + 2] = cc;
            if (edge_to_sites.contains(cc)) {
                edge_to_sites[cc].insert(i);
                edge_to_sites[cc].insert(i + 1);
                edge_to_sites[cc].insert(i + 2);
            } else
                edge_to_sites.insert({cc, {i, i + 1, i + 2}});
        }
        for (const auto& [vert, neighbors]: edge_to_sites) {
            for(auto i: neighbors) {
                if (m_delauney.halfedges[i] == delaunator::INVALID_INDEX) {
                    auto edge_index = m_voroni_edges.size();
                    RegionEdge edge{{vert, vert}};
                    m_voroni_edges.push_back(edge);

                    auto r0 = m_delauney.triangles[i];
                    auto r1 = m_delauney.triangles[i % 3 == 2 ? i - 2 : i + 1];
                    m_regions[r0].neighborhood[r1] = edge_index;
                    m_regions[r1].neighborhood[r0] = edge_index;
                } else if (i > m_delauney.halfedges[i]) {
                    const auto& next_vert = site_to_edge[m_delauney.halfedges[i]];
                    auto edge_index = m_voroni_edges.size();
                    RegionEdge edge{{vert, next_vert}};
                    m_voroni_edges.push_back(edge);

                    auto r0 = m_delauney.triangles[i];
                    auto r1 = m_delauney.triangles[m_delauney.halfedges[i]];
                    m_regions[r0].neighborhood[r1] = edge_index;
                    m_regions[r1].neighborhood[r0] = edge_index;
                }
            }
        }

        // set ocean to hull of map
        for(auto i = m_delauney.hull_start; m_delauney.hull_next[i] != m_delauney.hull_start; i = m_delauney.hull_next[i])
            m_oceans.insert(i);

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
                auto &region = m_regions[i];
                region.terrain.sea_level = level;
                for (const auto& [n, edge]: region.neighborhood) {
                    if (!added.contains(n)) {
                        next_layer.insert(n);
                        added.insert(n);
                    }
                }
            }
            layer = next_layer;
            level++;
        } while (!layer.empty());

        // randomly choose tectonic plate origins from faces
        added.clear();
        std::unordered_set<std::size_t> tectonic_regions;
        for (auto index = 0; index < m_regions.size(); index++) {
            m_regions[index].terrain.tectonic_plate = tectonic_regions.size();
            tectonic_regions.insert(index);
            added.insert(index);
            if(tectonic_regions.size() >= num_tectonic_plates)
                break;
        }

        // assign terrain regions to tectonic plates and make mountains along plate boundaries
        int mountain_size = 2;
        while (added.size() < m_regions.size()) {
            std::unordered_set<std::size_t> next;
            for (auto index : tectonic_regions) {
                const auto &region = m_regions[index];
                for (const auto& [n, edge]: region.neighborhood) {
                    if (!added.contains(n)) {
                        m_regions[n].terrain.tectonic_plate = region.terrain.tectonic_plate;
                        added.insert(n);
                        next.insert(n);
                    } else if (!m_oceans.contains(n) // region can't be mountain and ocean
                    && m_regions[n].terrain.tectonic_plate != region.terrain.tectonic_plate // mountains on tectonic plate boundaries
                    && m_regions[n].terrain.sea_level > mountain_size + 2) // no mountains on beach
                        m_mountains.insert(n);
                }
            }
            tectonic_regions = next;
        }

        // assign elevations based on distance from mountain or distance from ocean (whichever is closer)
        for (auto i = 0; i < m_regions.size(); i++) {
            auto& region = m_regions[i];
            if (m_oceans.contains(i)) {
                region.terrain.elevation = 0.0;
                region.terrain.humidity = 1.0;
            } else if (m_mountains.contains(i)) {
                region.terrain.elevation = 1.0;
                region.terrain.humidity = 0.0;
            } else {
                float f = region.terrain.sea_level/(level - 1.0f);
                region.terrain.elevation = f;
                region.terrain.humidity = 1.0f - f;
                continue;
                int path_length;
                find_nearest_mountain_face(i, path_length);
                if (path_length > mountain_size)
                    region.terrain.elevation = (0.2f * region.terrain.sea_level) / level;
                else {
                    auto d = (10.0 * (path_length - 1)) / mountain_size;
                    region.terrain.elevation = 0.8f / glm::log(d + glm::exp(1));
                }
            }
        }

        // wind
        auto mountain_shadow = 3*mountain_size; // # of layers after the mountain cell affected by mountain
        std::map<float, std::size_t> wind_projections;
        for (auto i = 0; i < m_regions.size(); i++) {
            break;
            const auto& region = m_regions[i];
            auto dx = glm::cos(wind_direction);
            auto dy = glm::sin(wind_direction);
            wind_projections[m_voroni_sites[i]*dx + m_voroni_sites[i+1]*dy] = i;
        }

        for (auto [wind_projection, i]: wind_projections) {
            break;
            auto& region = m_regions[i];
            if (m_oceans.contains(region.coord_index) || m_mountains.contains(region.coord_index))
                continue;
            int path_length;
            auto mountain_region = m_regions[find_nearest_mountain_face(i, path_length)];
            auto mountain_site = site_at(mountain_region.coord_index);
            auto dx = glm::cos(wind_direction);
            auto dy = glm::sin(wind_direction);
            auto mountain_projection = dx*mountain_site.x + dy*mountain_site.y;
            if (path_length > mountain_size) {
                if (path_length >= region.terrain.sea_level)
                    region.terrain.humidity = 1.f;
                else if (wind_projection < mountain_projection || path_length > mountain_shadow) {
                    float area_humidity = 0.0;
                    auto count = 0;
                    for (const auto& [n, edge]: region.neighborhood) {
                        auto neighbor = m_regions[n];
                        auto neighbor_site = site_at(neighbor.coord_index);
                        auto neighbor_projection = dx*neighbor_site.x + dy*neighbor_site.y;
                        if(neighbor_projection < wind_projection)
                            area_humidity += EVAPORATION*neighbor.terrain.humidity;
                    }
                    area_humidity = std::min(1.0f, area_humidity);
                    region.terrain.humidity = area_humidity;
                } else if (path_length < mountain_shadow) {
                    region.terrain.humidity = 0.2f + 0.2f*path_length/mountain_shadow;
                }
            } else
                region.terrain.humidity = 1.5f - region.terrain.elevation;
        }
    }

    entt::entity Terrain::register_mesh(entt::registry &registry) {
        m_bullet_mesh = std::make_shared<btTriangleMesh>();

        Mesh mesh;
        auto ocean_color = glm::vec3(0, 0, 1);
        auto desert_color = glm::vec3(251, 225, 182)/255.f;
        auto grassland_color = glm::vec3(83, 91, 41)/255.f;
        auto mountain_color = glm::vec3(.6, .6, .6);
        auto mountain_top_color = glm::vec3(1, 1, 1);
        for (std::size_t i = 0; i < m_delauney.triangles.size(); i += 3) {
            auto color = ocean_color;
            for(auto j = i; j < i + 3; j++) {
                auto index = m_delauney.triangles[j];
                const auto& region = m_regions[index];
                auto elevation = region.terrain.elevation;
                auto site = site_at(region.coord_index);

                if (!m_oceans.contains(index)) {
                    elevation *= MOUNTAIN_HEIGHT;
                    if(region.terrain.sea_level == 1 || region.terrain.humidity < .3)
                        color = desert_color;
                    else
                        color = grassland_color * region.terrain.humidity;
                    if (elevation > .6f * MOUNTAIN_HEIGHT)
                        color = mountain_color;
                    if (elevation >= .9f * MOUNTAIN_HEIGHT)
                        color = mountain_top_color;
                }

                mesh.vertices.emplace_back(site.x, elevation, site.y);
                mesh.colors.push_back(color);
                mesh.indices.push_back(j);
            }
            auto v1 = glm2bt(mesh.vertices[i]);
            auto v2 = glm2bt(mesh.vertices[i + 1]);
            auto v3 = glm2bt(mesh.vertices[i + 2]);
            m_bullet_mesh->addTriangle(v1, v2, v3);
        }
        // register physics data
        m_shape = std::make_shared<btBvhTriangleMeshShape>(m_bullet_mesh.get(), true);
        m_body = std::make_shared<btCollisionObject>();
        m_body->setCollisionShape(m_shape.get());

        // register mesh
        m_entity = registry.create();
        registry.emplace_or_replace<Mesh>(m_entity, mesh);
        registry.patch<InstanceList>(m_entity, [](auto &instance_list) {
            instance_list.set_instances(std::vector<glm::mat4>{glm::mat4(1)});
        });

        return m_entity;
    }

    std::size_t Terrain::find_nearest_mountain_face(std::size_t index, int& path_length) {
        assert(!m_mountains.empty());
        if(m_mountains.contains(index))
            return index;
        std::unordered_set<std::size_t> layer;
        std::unordered_set<std::size_t> searched{index};
        auto min_distance = std::numeric_limits<double>::infinity();
        const auto& start_region = m_regions[index];
        for(const auto& [n, edge]: start_region.neighborhood)
            layer.insert(n);
        auto found{index};
        int path_distance{0};
        while(found == index && searched.size() < m_regions.size() && layer.empty()) {
            path_distance++;
            std::unordered_set<std::size_t> next;
            for (auto i: layer) {
                if (m_mountains.contains(i)) {
                    auto distance = glm::distance(
                            site_at(m_regions[i].coord_index),
                            site_at(start_region.coord_index));
                    if (distance < min_distance) {
                        found = i;
                        min_distance = distance;
                    }
                } else {
                    for(const auto& [n, edge]: m_regions[i].neighborhood) {
                        if(!searched.contains(n))
                            next.insert(n);
                    }
                }
                searched.insert(i);
            }
            layer = next;
        }
        path_length = path_distance;
        return found;
    }

    glm::vec3 Terrain::get_mouse_terrain_collision_point(float x, float y, const RenderContext& context) const {
        if(m_bullet_mesh == nullptr || m_shape == nullptr)
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

    std::vector<glm::vec3> Terrain::get_mouse_terrain_collision_triangle(float x,
                                                                         float y,
                                                                         const RenderContext& context) const {
        if(m_bullet_mesh == nullptr || m_shape == nullptr)
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
            auto& tmesh = m_bullet_mesh->getIndexedMeshArray()[0];
            auto* gfxbase = (unsigned int*)(tmesh.m_triangleIndexBase+res.m_index*tmesh.m_triangleIndexStride);
            for (int j = 0; j < 3; j++) {
                int graphicsindex = tmesh.m_indexType==PHY_SHORT?((unsigned short*)gfxbase)[j]:gfxbase[j];
                auto* graphicsbase = (float*)(tmesh.m_vertexBase+graphicsindex*tmesh.m_vertexStride);
                btVector3 v(graphicsbase[0] * m_bullet_mesh->getScaling().getX(),
                            graphicsbase[1] * m_bullet_mesh->getScaling().getY(),
                            graphicsbase[2] * m_bullet_mesh->getScaling().getZ());
                verts.emplace_back(v.x(), v.y(), v.z());
            }

            return verts;
        }
        return {};
    }

    std::vector<double> Terrain::generate_sites(unsigned int num_regions, unsigned int width, unsigned int height) {
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

} // namespace mapgen
