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
    std::vector<Edge> edges;
};


#endif //MAPGEN_DIAGRAM_H
