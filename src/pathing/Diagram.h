//
// Created by Philip Smith on 4/11/2021.
//

#ifndef MAPGEN_DIAGRAM_H
#define MAPGEN_DIAGRAM_H

#include <glm/glm.hpp>
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

    // Convert from Voroni to Delauney or Delauney to Voroni
    Diagram dual() {
        Diagram diagram;

        for(const auto& edge: edges) {
            unsigned int index = diagram.faces.size() - 1;
            auto found{false};
            for(auto f: diagram.faces) {
                if(f.site == edge.first) {
                    f.neighbors.emplace_back(index + 1);
                    found = true;
                    break;
                }
            }
            if(!found)
                diagram.faces.emplace_back(edge.first, std::vector<unsigned int>{index + 1});
            index++;
            found = false;
            for(auto f: diagram.faces) {
                if(f.site == edge.first) {
                    f.neighbors.emplace_back(index - 1);
                    found = true;
                    break;
                }
            }
            if(!found)
                diagram.faces.emplace_back(edge.second, std::vector<unsigned int>{index - 1});
        }
        for(const auto& face: faces) {
            for(auto neighbor: face.neighbors) {
                diagram.edges.emplace_back(face.site, faces[neighbor].site);
            }
        }

        return diagram;
    }
};


#endif //MAPGEN_DIAGRAM_H
