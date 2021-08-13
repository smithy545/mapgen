//
// Created by Philip Smith on 4/17/2021.
//

#include <mapgen/Terrain.h>

#include <BulletCollision/NarrowPhaseCollision/btRaycastCallback.h>
#include <cassert>
#include <engine/InstanceList.h>
#include <glm/glm.hpp>
#include <glm/gtx/matrix_decompose.hpp>
#include <glm/gtx/projection.hpp>
#include <limits>
#include <nlohmann/json.hpp>
#include <random>
#include <stdexcept>
#include <utils/file_util.h>
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
                     : Terrain(generate_sites(num_regions, width, height),
                               num_tectonic_plates,
                               wind_direction) {}

    Terrain::Terrain(std::vector<double> site_coords, unsigned int num_tectonic_plates, float wind_direction)
    : m_delauney(site_coords) {
        assert(num_tectonic_plates > 0);


        // Make regions from site coords
        m_voroni_sites = std::move(site_coords);
        for(std::size_t i = 0; i < m_voroni_sites.size(); i += 2)
            m_regions.push_back(Region{i});


        // assign region relations/neighborhoods and generate voroni edges from delauney triangulation
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
                std::size_t edge_index{m_voroni_edges.size()};
                if (m_delauney.halfedges[i] == delaunator::INVALID_INDEX) {
                    RegionEdge edge{{vert, vert}};
                    m_voroni_edges.push_back(edge);

                    auto r0 = m_delauney.triangles[i];
                    auto r1 = m_delauney.triangles[i % 3 == 2 ? i - 2 : i + 1];
                    m_regions[r0].neighborhood[r1] = edge_index;
                    m_regions[r1].neighborhood[r0] = edge_index;
                } else if (i > m_delauney.halfedges[i]) {
                    const auto& next_vert = site_to_edge[m_delauney.halfedges[i]];
                    RegionEdge edge{{vert, next_vert}};
                    m_voroni_edges.push_back(edge);

                    auto r0 = m_delauney.triangles[i];
                    auto r1 = m_delauney.triangles[m_delauney.halfedges[i]];
                    m_regions[r0].neighborhood[r1] = edge_index;
                    m_regions[r1].neighborhood[r0] = edge_index;
                }
            }
        }
        export_to_file("recent.json");


        // set ocean to hull of map
        std::unordered_set<std::size_t> oceans;
        for(auto i = m_delauney.hull_start; m_delauney.hull_next[i] != m_delauney.hull_start; i = m_delauney.hull_next[i])
            oceans.insert(i);
        // assign ocean field (sea level = path length to map hull)
        std::size_t ocean_field_max = 0;
        m_ocean_field = std::make_shared<PathLengthField>(
            *this,
            std::vector<std::size_t>(oceans.begin(), oceans.end()),
            [&ocean_field_max, &oceans](
                    const PathLengthField& field,
                    const Terrain& terrain,
                    std::size_t index) {
                if(oceans.contains(index))
                    return PathStats{0, index, index};

                PathStats closest_ocean{INT_MAX, index, index};
                for(const auto& [n, e]: terrain.get_region(index).neighborhood) {
                    if(field.contains(n)) {
                        auto ocean = field[n];
                        if(ocean.length < closest_ocean.length) {
                            closest_ocean.length = ocean.length;
                            closest_ocean.origin = ocean.origin;
                        }
                    }
                }
                closest_ocean.length += 1;
                if(closest_ocean.length > ocean_field_max)
                    ocean_field_max = closest_ocean.length;
                return closest_ocean;
            }
        );


        // randomly choose tectonic plate origins from faces
        std::unordered_set<std::size_t> mountains;
        std::unordered_set<std::size_t> added, tectonic_regions;
        for (auto index = 0; index < m_regions.size(); index++) {
            m_regions[index].terrain.tectonic_plate = tectonic_regions.size();
            tectonic_regions.insert(index);
            added.insert(index);
            if(tectonic_regions.size() >= num_tectonic_plates)
                break;
        }
        // assign terrain regions to tectonic plates and make mountains along plate boundaries
        int mountain_size = 2;
        std::unordered_set<std::size_t> plate_regions;
        std::unordered_set<std::size_t> plate_edges;
        while (added.size() < m_regions.size()) {
            std::unordered_set<std::size_t> next;
            for (auto index : tectonic_regions) {
                const auto &region = m_regions[index];
                for (const auto& [n, edge]: region.neighborhood) {
                    auto& neighbor = m_regions[n];
                    if (!added.contains(n)) {
                        neighbor.terrain.tectonic_plate = region.terrain.tectonic_plate;
                        added.insert(n);
                        next.insert(n);
                    } else if (neighbor.terrain.tectonic_plate != region.terrain.tectonic_plate) {
                        plate_edges.insert(edge);
                        plate_regions.insert(index);
                        plate_regions.insert(n);
                        mountains.insert(n);
                    }
                }
            }
            tectonic_regions = next;
        }
        // field containing distance to nearest mountain
        auto mountain_field_max = 0;
        m_mountain_field = std::make_shared<PathLengthField>(
            *this,
            std::vector<std::size_t>(mountains.begin(), mountains.end()),
            [&mountain_field_max, &mountains](
                    const PathLengthField& field,
                    const Terrain& terrain,
                    std::size_t index) {
                if(mountains.contains(index))
                    return PathStats{0, index, index};

                PathStats closest_mountain{INT_MAX, index, index};
                for(const auto& [n, e]: terrain.get_region(index).neighborhood) {
                    if(field.contains(n)) {
                        auto mountain = field[n];
                        if(mountain.length < closest_mountain.length) {
                            closest_mountain.length = mountain.length;
                            closest_mountain.origin = mountain.origin;
                        }
                    }
                }
                closest_mountain.length += 1;
                if(closest_mountain.length > mountain_field_max)
                    mountain_field_max = closest_mountain.length;
                return closest_mountain;
            }
        );
        // assign elevations based on distance from mountain or distance from ocean (whichever is closer)
        for (auto i = 0; i < m_regions.size(); i++) {
            auto& region = m_regions[i];
            if (is_ocean(i)) {
                region.terrain.elevation = 0.0;
                region.terrain.precipitation = 1.0;
            } else if (is_mountain(i)) {
                region.terrain.elevation = 1.0;
                region.terrain.precipitation = 0.0;
            } else {
                auto path_to_closest_mountain = (*m_mountain_field)[i];
                auto path_to_closest_ocean = (*m_ocean_field)[i];
                if (path_to_closest_mountain.length > mountain_size)
                    region.terrain.elevation = (0.2 * path_to_closest_ocean.length) / ocean_field_max;
                else {
                    auto d = (10.0 * (path_to_closest_mountain.length - 1)) / mountain_size;
                    region.terrain.elevation = 0.8 / glm::log(d + glm::exp(1));
                }
                region.terrain.precipitation = 1.0f - region.terrain.elevation;
            }
        }


        // order regions by projection along wind vector
        glm::vec2 wind_vector(glm::cos(wind_direction), glm::sin(wind_direction));
        std::map<float, std::size_t> wind_projections;
        for (auto i = 0; i < m_regions.size(); i++)
            wind_projections[m_voroni_sites[i]*wind_vector.x + m_voroni_sites[i+1]*wind_vector.y] = i;

        std::vector<std::size_t> ordered_regions;
        for(const auto& [projection, index]: wind_projections)
            ordered_regions.push_back(index);

        m_wind_field = std::make_shared<FlatVectorField>(
            *this,
            ordered_regions,
            [&](
                    const TerrainField<glm::vec2>& field,
                    const Terrain& terrain,
                    std::size_t index) {
                if(terrain.is_ocean(index)) // oceans are wind sources
                    return wind_vector;
                if(terrain.is_mountain(index)) // mountains are wind sinks/breakers
                    return glm::vec2(0,0);

                glm::vec2 wind_direction(0,0);
                float num_neighbors = 0;
                for(const auto& [n, e]: terrain.get_region(index).neighborhood) {
                    if(field.contains(n)) {
                        auto dv = glm::normalize(terrain.site_at(index) - terrain.site_at(n));
                        wind_direction += glm::proj(dv, field[n]);
                        num_neighbors += 1;
                    }
                }
                wind_direction /= num_neighbors;
                wind_direction += wind_vector * m_wind_strength; // constant wind everywhere
                return wind_direction;
            }
        );


        // assign precipitation/humidity
        auto mountain_shadow = 3*mountain_size; // # of layers after the mountain cell affected by mountain
        for (auto [wind_projection, i]: wind_projections) {
            break;
            auto& region = m_regions[i];
            if (is_ocean(region.coord_index) || is_mountain(region.coord_index))
                continue;
            int path_length;
            auto mountain_region = m_regions[(*m_mountain_field)[i].origin];
            auto mountain_site = site_at(mountain_region.coord_index);
            auto dx = glm::cos(wind_direction);
            auto dy = glm::sin(wind_direction);
            auto mountain_projection = dx*mountain_site.x + dy*mountain_site.y;
            if (path_length > mountain_size) {
                if (path_length >= (*m_ocean_field)[i].length)
                    region.terrain.precipitation = 1.f;
                else if (wind_projection < mountain_projection || path_length > mountain_shadow) {
                    float area_humidity = 0.0;
                    auto count = 0;
                    for (const auto& [n, edge]: region.neighborhood) {
                        auto neighbor = m_regions[n];
                        auto neighbor_site = site_at(neighbor.coord_index);
                        auto neighbor_projection = dx*neighbor_site.x + dy*neighbor_site.y;
                        if(neighbor_projection < wind_projection)
                            area_humidity += EVAPORATION*neighbor.terrain.precipitation;
                    }
                    area_humidity = std::min(1.0f, area_humidity);
                    region.terrain.precipitation = area_humidity;
                } else if (path_length < mountain_shadow) {
                    region.terrain.precipitation = 0.2f + 0.2f*path_length/mountain_shadow;
                }
            } else
                region.terrain.precipitation = 1.5f - region.terrain.elevation;
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
            for(auto j = i; j < i + 3; j++) {
                auto t = m_delauney.triangles[j];
                const auto& terrain = m_regions[m_delauney.triangles[j]].terrain;
                auto elevation = terrain.elevation;
                auto color = ocean_color;
                if (!is_ocean(t)) {
                    elevation *= MOUNTAIN_HEIGHT;
                    if(terrain.precipitation < .3)
                        color = desert_color;
                    else
                        color = grassland_color * terrain.precipitation;
                    if (elevation > .6f * MOUNTAIN_HEIGHT)
                        color = mountain_color;
                    if (elevation >= .9f * MOUNTAIN_HEIGHT)
                        color = mountain_top_color;
                }

                mesh.vertices.emplace_back(m_voroni_sites[2*t], elevation, m_voroni_sites[2*t + 1]);
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

    void Terrain::export_to_file(const std::string& filename) {
        json data;
        data["num_sites"] = m_voroni_sites.size()/2;
        for(auto v: m_voroni_sites)
            data["site_coords"].push_back(v);
        utils::file::write_json_file("/Users/philipsmith/CLionProjects/civilwar/res/" + filename, data);
    }

    bool Terrain::is_mountain(std::size_t index) const {
        return m_mountain_field->contains(index) && (*m_mountain_field)[index].length == 0;
    }

    bool Terrain::is_ocean(std::size_t index) const {
        return m_ocean_field->contains(index) && (*m_ocean_field)[index].length == 0;
    }

    template<typename T>
    TerrainField<T>::TerrainField(const Terrain &terrain,
                                      const std::vector<std::size_t> &seeds,
                                      const std::function<T(const TerrainField &,
                                                            const Terrain &,
                                                            std::size_t)> &field,
                                      const std::function<bool(const TerrainField &,
                                                           const Terrain &,
                                                           std::size_t)> &predicate) {
        assert(!seeds.empty());
        std::unordered_set<std::size_t> current{seeds.begin(), seeds.end()};
        while(m_values.size() < terrain.get_num_regions() && !current.empty()) {
            std::unordered_set<std::size_t> next;
            for(auto i: current) {
                m_values[i] = field(*this, terrain, i);
                for(const auto& [n, e]: terrain.get_region(i).neighborhood) {
                    if(!current.contains(n) && !m_values.contains(n) && predicate(this, terrain, n))
                        next.insert(n);
                }
            }
            current = next;
        }
    }

    template <typename T>
    TerrainField<T>::TerrainField(const Terrain& terrain,
                                     const std::vector<std::size_t> &seeds,
                                     const std::function<T(const TerrainField&,
                                                           const Terrain& terrain,
                                                           std::size_t)> &field) {
        assert(!seeds.empty());
        std::unordered_set<std::size_t> current{seeds.begin(), seeds.end()};
        while(m_values.size() < terrain.get_num_regions()) {
            std::unordered_set<std::size_t> next;
            for(auto i: current) {
                m_values[i] = field(*this, terrain, i);
                for(const auto& [n, e]: terrain.get_region(i).neighborhood) {
                    if(!current.contains(n) && !m_values.contains(n))
                        next.insert(n);
                }
            }
            current = next;
        }
    }
} // namespace mapgen
