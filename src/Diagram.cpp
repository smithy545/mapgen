//
// Created by Philip Smith on 6/29/21.
//

#include <mapgen/Diagram.h>
#include <mapgen/DelaunatorAlgorithm.h>


namespace mapgen {
    void Diagram::add_edge(glm::vec2 origin, glm::vec2 destination, unsigned int face1, unsigned int face2) {
        auto index = m_edges.size();
        m_edges.emplace_back(origin, destination);
        m_faces[face1].neighboring_edges[face2] = index;
        m_faces[face2].neighboring_edges[face1] = index;
    }

    int Diagram::add_face(glm::vec2 site) {
        auto index = m_faces.size();
        m_faces.emplace_back(site);
        return index;
    }

    void Diagram::add_to_hull(unsigned int index) {
        m_hull.insert(index);
    }

    bool Diagram::on_hull(unsigned int index) {
        return m_hull.contains(index);
    }

    // Convert from Voroni to Delauney or Delauney to Voroni
    Diagram Diagram::dual() const {
        // TODO: translate hull information (hull consists of sites for voroni and edges for Delauney) to dual
        Diagram diagram;
        std::unordered_map<std::string, unsigned int> added_faces;
        for (auto edge: m_edges) {
            auto k1 = site_key(edge.first);
            auto k2 = site_key(edge.second);
            if (!added_faces.contains(k1))
                added_faces[k1] = diagram.add_face(edge.first);
            if (!added_faces.contains(k2))
                added_faces[k2] = diagram.add_face(edge.second);
        }
        for (int i = 0; i < m_faces.size(); i++) {
            auto face = m_faces[i];
            for (auto[neighbor, edge]: face.neighboring_edges) {
                auto e = m_edges[edge];
                auto k1 = site_key(e.first);
                auto k2 = site_key(e.second);
                if (i > neighbor)
                    diagram.add_edge(face.site, m_faces[neighbor].site, added_faces[k1], added_faces[k2]);
            }
        }

        return diagram;
    }

    // apply lloyd relaxation
    Diagram Diagram::relax() const {
        Diagram diagram;
        std::vector<double> new_sites;
        for(auto face: m_faces) {
            glm::vec2 avg{0,0};
            for(auto [n, e]: face.neighboring_edges)
                avg += m_edges[e].first;
            avg /= face.neighboring_edges.size();
            new_sites.push_back(avg.x);
            new_sites.push_back(avg.y);
        }
        return DelaunatorAlgorithm::construct_voroni(new_sites);
    }

    std::string Diagram::edge_key(glm::vec2 origin, glm::vec2 destination) {
        return fmt::format("{}, {} | {}, {}", origin.x, origin.y, destination.x, destination.y);
    }

    std::string Diagram::site_key(glm::vec2 site, bool plain) {
        if (plain)
            return fmt::format("{},{}", (int) site.x, (int) site.y);
        return fmt::format("glm::vec2({},{})", (int) site.x, (int) site.y);
    }
} // namespace mapgen