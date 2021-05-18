//
// Created by Philip Smith on 4/11/2021.
//

#ifndef MAPGEN_DIAGRAM_H
#define MAPGEN_DIAGRAM_H

#include <fmt/format.h>
#include <glm/glm.hpp>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <utils/macros.h>
#include <utility>
#include <vector>

namespace mapgen {
    struct Diagram {
        typedef std::pair<glm::vec2, glm::vec2> Edge;

        struct Face {
            typedef std::shared_ptr<Face> Ptr;

            Face(glm::vec2 site) : site(site), neighboring_edges() {}

            glm::vec2 site;
            std::unordered_map<unsigned int, unsigned int> neighboring_edges;
        };

        void add_edge(glm::vec2 origin, glm::vec2 destination, unsigned int face1, unsigned int face2) {
            int index = m_edges.size();
            m_edges.emplace_back(origin, destination);
            m_faces[face1].neighboring_edges[face2] = index;
            m_faces[face2].neighboring_edges[face1] = index;
        }

        int add_face(glm::vec2 site) {
            int index = m_faces.size();
            m_faces.emplace_back(site);
            return index;
        }

        void add_to_hull(unsigned int index) {
            m_hull.insert(index);
        }

        bool on_hull(unsigned int index) {
            return m_hull.contains(index);
        }

        // Convert from Voroni to Delauney or Delauney to Voroni
        Diagram dual() const {
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

        static std::string edge_key(glm::vec2 origin, glm::vec2 destination) {
            return fmt::format("{}, {} | {}, {}", origin.x, origin.y, destination.x, destination.y);
        }

        static std::string site_key(glm::vec2 site, bool plain = true) {
            if (plain)
                return fmt::format("{},{}", (int) site.x, (int) site.y);
            return fmt::format("glm::vec2({},{})", (int) site.x, (int) site.y);
        }

    REFVAR_GET(std::vector<Edge>, edges, public);
    REFVAR_GET(std::vector<Face>, faces, public);
    REFVAR_GET(std::unordered_set<unsigned int>, hull, public);
    };
} // namespace mapgen

#endif //MAPGEN_DIAGRAM_H
