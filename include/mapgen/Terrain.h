//
// Created by Philip Smith on 4/17/2021.
//

#ifndef MAPGEN_TERRAIN_H
#define MAPGEN_TERRAIN_H

#include <btBulletCollisionCommon.h>
#include <engine/mesh.h>
#include <engine/RenderContext.h>
#include <entt/entt.hpp>
#include <functional>
#include <mapgen/Diagram.h>
#include <string>
#include <unordered_set>
#include <utils/macros.h>
#include <vector>


using namespace engine;

namespace mapgen {
    struct TerrainRegion {
        const Diagram::Face& face;
        float elevation{0.0};
        float moisture{0.0};
        int plate_id{-1};
        int sea_level{0};
    };

    class Terrain {
    public:
        PTR(Terrain);

        Terrain(unsigned int num_sites,
                int width,
                int height,
                int num_tectonic_plates = 1,
                glm::vec2 wind_direction = glm::vec2(1,0));

        void register_mesh(entt::registry &registry);

        entt::entity register_wireframe_mesh(entt::registry &registry);

        entt::entity register_voroni_mesh(entt::registry &registry);

        [[nodiscard]]
        glm::vec3 get_mouse_terrain_collision(float x, float y, const RenderContext& context) const;

        [[nodiscard]]
        std::vector<glm::vec3> get_mouse_terrain_colliding_triangle(float x, float y, const RenderContext& context) const;

        [[nodiscard]]
        Mesh get_region_mesh(unsigned int index);

    private:
        std::unordered_set<unsigned int> m_mountains;
        std::unordered_set<unsigned int> m_oceans;
        std::unordered_set<unsigned int> m_rivers;
        std::vector<TerrainRegion> m_regions;
        std::shared_ptr<btBvhTriangleMeshShape> m_shape{nullptr};
        std::shared_ptr<btCollisionObject> m_body{nullptr};
        std::shared_ptr<btTriangleMesh> m_mesh{nullptr};

        void assign_ocean_from_hull(int depth);

        void apply_callback_breadth_first(const std::function<void(unsigned int)>& callback,
                                        unsigned int start,
                                        int max_depth,
                                        const std::function<bool(unsigned int)>& predicate =
                                                [](unsigned int){return true;},
                                        std::unordered_set<unsigned int> visited = {});

        void apply_callback_breadth_first_recursively(const std::function<void(unsigned int)>& callback,
                                                      unsigned int start,
                                                      int max_depth,
                                                      const std::function<bool(unsigned int)>& predicate =
                                                              [](unsigned int){return true;},
                                                      std::unordered_set<unsigned int> visited = {});

        unsigned int find_nearest_mountain_face(unsigned int index, int& path_length);

    VAR_GET(entt::entity, entity, public){entt::null};
    VAR_GET(Diagram, base, public);
    VAR_GET(Diagram, dual, public);
    };
} // namespace mapgen

#endif //MAPGEN_TERRAIN_H
