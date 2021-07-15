//
// Created by Philip Smith on 4/19/2021.
//

#ifndef MAPGEN_DELAUNATORALGORITHM_H
#define MAPGEN_DELAUNATORALGORITHM_H

#include <delaunator-header-only.hpp>
#include <stdexcept>
#include <vector>

#include "Diagram.h"


namespace mapgen {
    class DelaunatorAlgorithm {
    public:
        static Diagram construct_voroni_diagram(const std::vector<double> &coords);

        static Diagram construct_clamped_voroni_diagram(const std::vector<double> &coords, double x, double y, double width, double height);

        static Diagram construct_delauney_diagram(const std::vector<double> &coords);
    private:
        static bool collides(glm::vec2 p, double min_x, double min_y, double max_x, double max_y) {
            return p.x >= min_x && p.x < max_x && p.y >= min_y && p.y < max_y;
        }

        static double cross(glm::vec2 v, glm::vec2 w) {
            return v.x*w.y - v.y*w.x;
        }

        static glm::vec2 compute_intersection(glm::vec2 p1, glm::vec2 p2, glm::vec2 q1, glm::vec2 q2) {
            // Reference: https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect/565282#565282
            auto r = p2 - p1;
            auto s = q2 - q1;
            auto d = q1 - p1;
            auto c1 = cross(r, s);
            auto c2 = cross(d, r);
            auto c3 = cross(d, s);
            if(c1 == 0) {
                if(c2 == 0) {
                    float rr = glm::dot(r, r);
                    float t0 = glm::dot(d, r) / rr;
                    float t1 = t0 + glm::dot(s, r)  / rr;
                    if((t0 >= 0 && t0 <= 1) || (t1 >= 0 && t1 <= 1))
                        return p1 + t0*r;
                }
            } else {
                float u = c2 / c1;
                float t = c3 / c1;
                if(t >= 0 && t <=1 && u >= 0 && u <= 1)
                    return p1 + t*r;
            }
            throw std::runtime_error("Lines do not intersect");
        }
    };
} // namespace mapgen

#endif //MAPGEN_DELAUNATORALGORITHM_H
