//
// Created by Philip Smith on 4/17/2021.
//

#ifndef MAPGEN_TERRAIN_H
#define MAPGEN_TERRAIN_H

#include <btBulletCollisionCommon.h>
#include <delaunator.hpp>
#include <engine/mesh.h>
#include <engine/RenderContext.h>
#include <entt/entt.hpp>
#include <fmt/format.h>
#include <functional>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utils/macros.h>
#include <utils/math_util.h>
#include <vector>


using namespace engine;

namespace mapgen {
    struct RegionTerrain {
        float elevation{0.0};
        float humidity{0.0};
        int tectonic_plate{-1};
        int sea_level{0};
    };

    struct RegionEdge {
        glm::vec2 vertices[2];

        inline bool is_full_edge() {
            return vertices[0] != vertices[1];
        }
    };

    struct Region {
        std::size_t coord_index;
        entt::entity id{entt::null};
        std::map<std::size_t, std::size_t> neighborhood{};
        RegionTerrain terrain{};
    };

    class Terrain {
    public:
        PTR(Terrain);

        Terrain(unsigned int num_regions,
                unsigned int width,
                unsigned int height,
                unsigned int num_tectonic_plates = 1,
                float wind_direction = .0f);

        explicit Terrain(std::vector<double> site_coords,
                unsigned int num_tectonic_plates = 1,
                float wind_direction = .0f);

        entt::entity register_mesh(entt::registry &registry);

        [[nodiscard]]
        glm::vec3 get_mouse_terrain_collision_point(float x, float y, const RenderContext& context) const;

        [[nodiscard]]
        std::vector<glm::vec3> get_mouse_terrain_collision_triangle(float x, float y, const RenderContext& context) const;

    private:
        const float MOUNTAIN_HEIGHT{200};
        const float EVAPORATION{0.3f};
        float m_wind_strength{0.1f}; // determines amount of evaporated moisture passing between neighboring regions
        float m_evaporation{0.1f};   // determines amount of moisture evaporating per region
        std::unordered_set<std::size_t> m_mountains;
        std::unordered_set<std::size_t> m_oceans;
        std::unordered_set<std::size_t> m_rivers;
        std::shared_ptr<btBvhTriangleMeshShape> m_shape{nullptr};
        std::shared_ptr<btCollisionObject> m_body{nullptr};
        std::shared_ptr<btTriangleMesh> m_bullet_mesh{nullptr};

        delaunator::Delaunator m_delauney;
        std::vector<RegionEdge> m_voroni_edges;
        std::vector<double> m_voroni_sites;

        std::size_t find_nearest_mountain_face(std::size_t index, int& path_length);

        inline glm::vec2 site_at(std::size_t index) {
            return {m_voroni_sites[index], m_voroni_sites[index + 1]};
        }

        static std::vector<double> generate_sites(unsigned int num_regions, unsigned int width, unsigned int height);

    VAR_GET(entt::entity, entity, public){entt::null};
    VAR_GET(std::vector<Region>, regions, public);
    };
} // namespace mapgen

#endif //MAPGEN_TERRAIN_H
