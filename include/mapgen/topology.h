//
// Created by Philip Smith on 7/26/21.
//

#ifndef MAPGEN_TOPOLOGY_H
#define MAPGEN_TOPOLOGY_H

#include <delaunator-header-only.hpp>
#include <functional>
#include <glm/glm.hpp>
#include <map>
#include <unordered_set>
#include <utility>
#include <utils/macros.h>
#include <vector>


namespace mapgen {
    class TerrainTopology;

    class RegionTopology {
    public:
        RegionTopology(std::size_t index, const TerrainTopology& parent);

        void add_neighbor(std::size_t neighbor);

        [[nodiscard]] glm::vec2 get_site() const;

        [[nodiscard]] std::size_t get_index() const;

        [[nodiscard]] std::vector<std::size_t> get_neighborhood() const;
    private:
        std::size_t m_index;
        const TerrainTopology& m_global;
        std::unordered_set<std::size_t> m_neighborhood{};
    };

    class TerrainTopology {
    public:
        explicit TerrainTopology(std::vector<double> coords);

        [[nodiscard]] std::vector<std::size_t> get_hull_regions() const;

        [[nodiscard]] std::size_t internal_to_region_index(std::size_t index) const;

        RegionTopology& get_region_topology(std::size_t index);
    private:
        std::map<std::size_t, std::size_t> m_region_map;

    REFVAR_GET(std::vector<double>, coords, public);
    REFVAR_GET(std::vector<RegionTopology>, regions, public);
    REFVAR_GET(delaunator::Delaunator, geometry, public);
    REFVAR_GET(std::vector<glm::vec2>, vertices, public);
    };
} // namespace mapgen


#endif //MAPGEN_TOPOLOGY_H
