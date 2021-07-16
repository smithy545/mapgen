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
            PTR(Face);

            explicit Face(glm::vec2 site) : site(site), neighboring_edges() {}

            glm::vec2 site;
            // edge indices in a map indexed by the corresponding neighbors face index
            std::unordered_map<unsigned int, unsigned int> neighboring_edges;
        };

        void add_edge(glm::vec2 origin, glm::vec2 destination, unsigned int face1, unsigned int face2);

        int add_face(glm::vec2 site);

        void add_to_hull(unsigned int index);

        bool on_hull(unsigned int index);

        void cull_edge(unsigned int index);

        void modify_edge(unsigned int index, glm::vec2 first, glm::vec2 second);

        Diagram dual() const;

        Diagram relax() const;

        Diagram relax(float x, float y, float width, float height) const;

        static std::string edge_key(glm::vec2 origin, glm::vec2 destination);

        static std::string site_key(glm::vec2 site, bool plain = true);

    REFVAR_GET(std::vector<Edge>, edges, public);
    REFVAR_GET(std::vector<Face>, faces, public);
    REFVAR_GET(std::unordered_set<unsigned int>, hull, public);
    };
} // namespace mapgen

#endif //MAPGEN_DIAGRAM_H
