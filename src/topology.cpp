//
// Created by Philip Smith on 7/26/21.
//

#include <mapgen/topology.h>

#include <utils/math_util.h>

#include <utility>


namespace mapgen {
    TerrainTopology::TerrainTopology(std::vector<double> coords) : m_coords(std::move(coords)), m_geometry(m_coords) {
        for(std::size_t e = 0; e < m_geometry.triangles.size(); e++) {
            auto i = m_geometry.triangles[e];
            if(!m_region_map.contains(i)) {
                m_region_map.insert({i, m_regions.size()});
                m_regions.emplace_back(i, *this);
            }
            auto& region = m_regions[m_region_map[i]];
            if(m_geometry.halfedges[e] != delaunator::INVALID_INDEX)
                region.add_neighbor(m_geometry.triangles[m_geometry.halfedges[e]]);
            if(e % 3 == 0) {
                auto i0 = 2*i;
                auto i1 = 2*m_geometry.triangles[e + 1];
                auto i2 = 2*m_geometry.triangles[e + 2];
                m_vertices.push_back(utils::math::compute_triangle_circumcenter(
                        glm::vec2(m_coords[i0], m_coords[i0 + 1]),
                        glm::vec2(m_coords[i1], m_coords[i1 + 1]),
                        glm::vec2(m_coords[i2], m_coords[i2 + 1])
                ));
            }
        }
    }

    std::vector<std::size_t> TerrainTopology::get_hull_regions() const {
        std::vector<std::size_t> regions;
        auto e = m_geometry.hull_start;
        do {
            regions.push_back(m_region_map.at(e));
            e = m_geometry.hull_next[e];
        } while(e != m_geometry.hull_start);
        return regions;
    }

    std::size_t TerrainTopology::internal_to_region_index(std::size_t index) const {
        return m_region_map.at(index);
    }

    RegionTopology& TerrainTopology::get_region_topology(std::size_t index) {
        return m_regions.at(index);
    }

    RegionTopology::RegionTopology(std::size_t index, const TerrainTopology &parent)
    :  m_index(index), m_global(parent) {}

    glm::vec2 RegionTopology::get_site() const {
        const auto& coords = m_global.get_coords();
        return glm::vec2(coords[2*m_index], coords[2*m_index + 1]);
    }

    std::size_t RegionTopology::get_index() const {
        return m_global.internal_to_region_index(m_index);
    }

    std::vector<std::size_t> RegionTopology::get_neighborhood() const {
        std::vector<std::size_t> regions;
        for(auto index: m_neighborhood)
            regions.push_back(m_global.internal_to_region_index(index));
        return regions;
    }

    void RegionTopology::add_neighbor(std::size_t neighbor) {
        m_neighborhood.insert(neighbor);
    }
} // namespace mapgen
