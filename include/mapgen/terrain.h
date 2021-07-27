//
// Created by Philip Smith on 4/17/2021.
//

#ifndef MAPGEN_TERRAIN_H
#define MAPGEN_TERRAIN_H

#include <btBulletCollisionCommon.h>
#include <engine/mesh.h>
#include <engine/RenderContext.h>
#include <entt/entt.hpp>
#include <fmt/format.h>
#include <functional>
#include <string>
#include <unordered_set>
#include <utils/macros.h>
#include <utils/math_util.h>
#include <vector>

#include "topology.h"


using namespace engine;

namespace mapgen {
    struct RegionTerrain {
        explicit RegionTerrain(const RegionTopology& neighbor_topology);

        [[nodiscard]] std::size_t get_index() const;

    VAR_GET(const RegionTopology&, topology, public);
    VAR(float, elevation, public, public){0.0};
    VAR(float, moisture, public, public){0.0};
    VAR(int, plate, public, public){-1};
    VAR(int, sea_level, public, public){0};
    };

    class Terrain {
    public:
        PTR(Terrain);

        Terrain(unsigned int num_regions,
                unsigned int width,
                unsigned int height,
                unsigned int num_tectonic_plates = 1,
                float wind_angle = .0f);

        explicit Terrain(std::vector<double> coords,
                         unsigned int num_tectonic_plates = 1,
                         float wind_angle = .0f);

        entt::entity register_mesh(entt::registry &registry);

        entt::entity register_wireframe_mesh(entt::registry &registry);

        [[nodiscard]]
        glm::vec3 get_mouse_terrain_collision(float x, float y, const RenderContext& context) const;

        [[nodiscard]]
        std::vector<glm::vec3> get_mouse_terrain_colliding_triangle(float x, float y, const RenderContext& context) const;

    private:
        const float MOUNTAIN_HEIGHT{200};
        float m_wind_strength{0.1}; // determines amount of evaporated moisture passing between neighboring regions
        float m_evaporation{0.1};   // determines amount of moisture evaporating per region
        std::unordered_set<std::size_t> m_mountains;
        std::unordered_set<std::size_t> m_oceans;
        std::unordered_set<std::size_t> m_rivers;
        std::vector<RegionTerrain> m_regions;
        std::shared_ptr<btBvhTriangleMeshShape> m_shape{nullptr};
        std::shared_ptr<btCollisionObject> m_body{nullptr};
        std::shared_ptr<btTriangleMesh> m_mesh{nullptr};

        void assign_ocean_from_hull(int depth);

        void apply_callback_breadth_first(const std::function<void(RegionTerrain&, std::size_t)>& callback,
                                          RegionTerrain& start,
                                          int max_depth,
                                          const std::function<bool(RegionTerrain&)>& predicate =
                                                [](RegionTerrain&){return true;},
                                          std::unordered_set<std::size_t> visited = {});

        void apply_callback_sweep_line(const std::function<void(RegionTerrain&, float)>& callback,
                                       float direction_in_radians,
                                       const std::function<bool(RegionTerrain&, float)>& predicate =
                                               [](RegionTerrain&, float){return true;});

        std::size_t find_nearest_mountain_face(std::size_t index, int& path_length);

        static std::vector<double> generate_sites(unsigned int num_regions, unsigned int width, unsigned int height);

    VAR_GET(entt::entity, entity, public){entt::null};
    VAR_GET(TerrainTopology, topology, public);
    };
} // namespace mapgen

#endif //MAPGEN_TERRAIN_H
