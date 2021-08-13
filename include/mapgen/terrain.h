//
// Created by Philip Smith on 4/17/2021.
//

#ifndef MAPGEN_TERRAIN_H
#define MAPGEN_TERRAIN_H

#include <btBulletCollisionCommon.h>
#include <delaunator.hpp>
#include <engine/mesh/Mesh.h>
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
        float elevation{0};     // meters
        float precipitation{0}; // centimeters
        float temperature{0};   // centigrade
        int tectonic_plate{-1};
    };

    struct RegionEdge {
        glm::vec2 vertices[2];

        [[nodiscard]] inline bool is_full_edge() const {
            return vertices[0] != vertices[1];
        }
    };

    struct Region {
        std::size_t coord_index;
        entt::entity id{entt::null};
        std::map<std::size_t, std::size_t> neighborhood{};
        RegionTerrain terrain{};
    };

    class Terrain;

    template <typename ValueType>
    class TerrainField {
    public:
        PTR(TerrainField);

        typedef std::map<std::size_t, ValueType> ValueMap;

        TerrainField(const Terrain& terrain,
                     const std::vector<std::size_t>& seeds,
                     const std::function<ValueType(const TerrainField&,
                                                   const Terrain&,
                                                   std::size_t)>& field,
                                                   const std::function<bool(const TerrainField&,
                                                           const Terrain&,
                                                           std::size_t)>& predicate);

        TerrainField(const Terrain& terrain,
                     const std::vector<std::size_t>& seeds,
                     const std::function<ValueType(const TerrainField&, const Terrain&, std::size_t)>& field);

        inline const ValueType& operator[](std::size_t index) const {
            return m_values.at(index);
        }

        [[nodiscard]] inline bool contains(std::size_t index) const {
            return m_values.contains(index);
        }

        REFVAR_GET(ValueMap, values, public){};
    };


    struct PathStats {
        int length;
        std::size_t origin;
        std::size_t destination;
    };
    typedef TerrainField<PathStats> PathLengthField;
    typedef TerrainField<glm::vec2> FlatVectorField;


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

        void export_to_file(const std::string& filename);

        [[nodiscard]]
        glm::vec3 get_mouse_terrain_collision_point(float x, float y, const RenderContext& context) const;

        [[nodiscard]]
        std::vector<glm::vec3> get_mouse_terrain_collision_triangle(float x, float y, const RenderContext& context) const;

        [[nodiscard]] bool is_mountain(std::size_t index) const;

        [[nodiscard]] bool is_ocean(std::size_t index) const;

        [[nodiscard]] inline const Region& get_region(std::size_t index) const {
            assert(index < m_regions.size());
            return m_regions[index];
        }

        [[nodiscard]] inline std::size_t get_num_regions() const {
            return m_regions.size();
        }

        [[nodiscard]] inline glm::vec2 site_at(std::size_t index) const {
            return {m_voroni_sites[index], m_voroni_sites[index + 1]};
        }

    private:
        const float MOUNTAIN_HEIGHT{200};
        const float EVAPORATION{0.3f};
        float m_wind_strength{0.1f}; // determines amount of evaporated moisture passing between neighboring regions
        float m_evaporation{0.1f};   // determines amount of moisture evaporating per region
        PathLengthField::Ptr m_mountain_field{nullptr};
        PathLengthField::Ptr m_ocean_field{nullptr};
        FlatVectorField::Ptr m_wind_field{nullptr};
        std::shared_ptr<btBvhTriangleMeshShape> m_shape{nullptr};
        std::shared_ptr<btCollisionObject> m_body{nullptr};
        std::shared_ptr<btTriangleMesh> m_bullet_mesh{nullptr};

        delaunator::Delaunator m_delauney;
        std::vector<RegionEdge> m_voroni_edges;
        std::vector<double> m_voroni_sites;

        static std::vector<double> generate_sites(unsigned int num_regions, unsigned int width, unsigned int height);

    VAR_GET(entt::entity, entity, public){entt::null};
    VAR_GET(std::vector<Region>, regions, public);
    };
} // namespace mapgen

#endif //MAPGEN_TERRAIN_H
