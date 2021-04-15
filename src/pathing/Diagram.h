//
// Created by Philip Smith on 4/11/2021.
//

#ifndef MAPGEN_DIAGRAM_H
#define MAPGEN_DIAGRAM_H

#include <fmt/format.h>
#include <glm/glm.hpp>
#include <unordered_set>
#include <utility>
#include <vector>


struct Diagram {
    typedef std::pair<glm::vec2, glm::vec2> Edge;

    struct Face {
        explicit Face(glm::vec2 site) : site(site), neighbors() {}
        Face(glm::vec2 site, std::vector<unsigned int> neighbors) : site(site), neighbors(std::move(neighbors)){}

        glm::vec2 site;
        std::vector<unsigned int> neighbors;
    };

    std::vector<Edge> edges;
    std::vector<Face> faces;

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
        for(const auto& edge: edges) {
            unsigned int index = diagram.faces.size();
            diagram.faces.emplace_back(edge.first, std::vector<unsigned int>{index + 1});
            diagram.faces.emplace_back(edge.second, std::vector<unsigned int>{index});
        }
        std::unordered_set<unsigned int> added_faces;
        for(int index = 0; index < faces.size(); index++) {
            auto face = faces[index];
            for(auto neighbor: face.neighbors) {
                if(!added_faces.contains(neighbor)) {
                    diagram.edges.emplace_back(face.site, faces[neighbor].site);
                }
            }
            added_faces.insert(index);
        }

        return diagram;
    }
};


#endif //MAPGEN_DIAGRAM_H
