//
// Created by Philip Smith on 4/11/2021.
//

#ifndef MAPGEN_DIAGRAM_H
#define MAPGEN_DIAGRAM_H

#include <fmt/format.h>
#include <glm/glm.hpp>
#include <memory>
#include <unordered_set>
#include <utils/macros.h>
#include <utility>
#include <vector>


class Diagram {
public:
    typedef std::pair<glm::vec2, glm::vec2> Edge;

    struct Face {
        typedef std::shared_ptr<Face> Ptr;
        explicit Face(glm::vec2 site) : site(site), neighbors() {}
        Face(glm::vec2 site, std::unordered_set<unsigned int> neighbors) : site(site), neighbors(std::move(neighbors)){}

        glm::vec2 site;
        std::unordered_set<unsigned int> neighbors;
    };

    void add_edge(glm::vec2 origin, glm::vec2 destination)  {
        m_edges.emplace_back(origin, destination);
    }

    void add_face(glm::vec2 site) {
        m_faces.push_back(Face{site});
    }

    void add_face(glm::vec2 site, std::unordered_set<unsigned int> neighbors) {
        m_faces.push_back(Face{site, neighbors});
    }

    void add_neighbor(unsigned int face, unsigned int neighbor) {
        m_faces[face].neighbors.insert(neighbor);
    }

    static std::string edge_key(glm::vec2 origin, glm::vec2 destination) {
        return fmt::format("{}, {} | {}, {}", origin.x, origin.y, destination.x, destination.y);
    }

    static std::string site_key(glm::vec2 site, bool plain = true) {
        if(plain)
            return fmt::format("{},{}", (int)site.x, (int)site.y);
        return fmt::format("glm::vec2({},{})", (int)site.x, (int)site.y);
    }

    // Convert from Voroni to Delauney or Delauney to Voroni
    Diagram dual() const {
        Diagram diagram;
        for(const auto& edge: m_edges) {
            unsigned int index = diagram.get_faces().size();
            diagram.add_face(edge.first, std::unordered_set<unsigned int>{index + 1});
            diagram.add_face(edge.second, std::unordered_set<unsigned int>{index});
        }
        std::unordered_set<unsigned int> added_faces;
        for(int index = 0; index < m_faces.size(); index++) {
            for(auto neighbor: m_faces[index].neighbors) {
                if(!added_faces.contains(neighbor))
                    diagram.add_edge(m_faces[index].site, m_faces[neighbor].site);
            }
            added_faces.insert(index);
        }

        return diagram;
    }

    REFVAR_GET(std::vector<Edge>, edges, public);
    REFVAR_GET(std::vector<Face>, faces, public);
};

#endif //MAPGEN_DIAGRAM_H
