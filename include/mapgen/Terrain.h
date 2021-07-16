//
// Created by Philip Smith on 4/17/2021.
//

#ifndef MAPGEN_TERRAIN_H
#define MAPGEN_TERRAIN_H

#include <btBulletCollisionCommon.h>
#include <engine/Mesh.h>
#include <engine/RenderContext.h>
#include <entt/entt.hpp>
#include <map>
#include <mapgen/Diagram.h>
#include <string>
#include <unordered_set>
#include <utils/macros.h>


using namespace engine;

namespace mapgen {
    struct TerrainRegion {
        unsigned int face;
        double elevation{0.0};
        unsigned int plate;
    };

    class Terrain {
    public:
        PTR(Terrain);

        Terrain(unsigned int num_sites, int width, int height, int num_tectonic_plates = 1, bool centered = false);

        void register_terrain_mesh(entt::registry &registry);

        void register_voroni_debug_mesh(entt::registry &registry);

        [[nodiscard]]
        glm::vec3 get_mouse_terrain_collision(float x, float y, const RenderContext& context) const;

        [[nodiscard]]
        std::vector<glm::vec3> get_mouse_terrain_colliding_triangle(float x, float y, const RenderContext& context) const;

        [[nodiscard]]
        Mesh get_face_mesh(unsigned int index, glm::vec3 color);

    private:
        std::unordered_set<unsigned int> m_mountains;
        std::unordered_set<unsigned int> m_oceans;
        std::unordered_set<unsigned int> m_rivers;
        std::unordered_set<unsigned int> m_tetonic_plates;
        std::map<unsigned int, TerrainRegion> m_regions;
        std::shared_ptr<btBvhTriangleMeshShape> m_shape{nullptr};
        std::shared_ptr<btCollisionObject> m_body{nullptr};
        std::shared_ptr<btTriangleMesh> m_mesh{nullptr};

        void assign_ocean(unsigned int start, int neighbor_depth);

        unsigned int find_nearest_mountain_face(unsigned int index);

        // TODO: add parameter to recursive search to stop searching after reaching neighbors 'X' hops away from index
        unsigned int find_nearest_mountain_face_recursive(unsigned int index);

        unsigned int find_mountain_kernel(unsigned int index, const std::unordered_set<unsigned int> &to_search,
                                          std::unordered_set<unsigned int> searched);

    VAR_GET(entt::entity, entity, public){entt::null};
    VAR_GET(Diagram, base, public);
    VAR_GET(Diagram, dual, public);
    };
} // namespace mapgen

#endif //MAPGEN_TERRAIN_H
