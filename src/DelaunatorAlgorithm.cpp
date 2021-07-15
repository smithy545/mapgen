//
// Created by Philip Smith on 4/19/2021.
//

#include <mapgen/DelaunatorAlgorithm.h>
#include <utils/math_util.h>


namespace mapgen {
    Diagram DelaunatorAlgorithm::construct_voroni_diagram(const std::vector<double> &coords) {
        Diagram voroni;
        delaunator::Delaunator d(coords);
        std::unordered_map<std::string, unsigned int> added_faces;
        std::unordered_map<unsigned int, glm::vec2> hull;
        for (std::size_t i = 0; i < d.triangles.size(); i++) {
            auto j = d.halfedges[i];
            auto i0 = i - (i % 3);
            auto c1 = utils::math::compute_triangle_circumcenter(
                    glm::vec2(d.coords[2 * d.triangles[i0]], d.coords[2 * d.triangles[i0] + 1]),
                    glm::vec2(d.coords[2 * d.triangles[i0 + 1]], d.coords[2 * d.triangles[i0 + 1] + 1]),
                    glm::vec2(d.coords[2 * d.triangles[i0 + 2]], d.coords[2 * d.triangles[i0 + 2] + 1])
            );
            glm::vec2 site(d.coords[2 * d.triangles[i]], d.coords[2 * d.triangles[i] + 1]);
            auto k = Diagram::site_key(site);
            if (!added_faces.contains(k))
                added_faces[k] = voroni.add_face(site);
            if (j != delaunator::INVALID_INDEX) {
                if (i > j) {
                    auto j0 = j - (j % 3);
                    auto p2 = glm::vec2(d.coords[2 * d.triangles[j]], d.coords[2 * d.triangles[j] + 1]);
                    auto c2 = utils::math::compute_triangle_circumcenter(
                            glm::vec2(d.coords[2 * d.triangles[j0]], d.coords[2 * d.triangles[j0] + 1]),
                            glm::vec2(d.coords[2 * d.triangles[j0 + 1]], d.coords[2 * d.triangles[j0 + 1] + 1]),
                            glm::vec2(d.coords[2 * d.triangles[j0 + 2]], d.coords[2 * d.triangles[j0 + 2] + 1])
                    );
                    voroni.add_edge(c1, c2, added_faces[k], added_faces[Diagram::site_key(p2)]);
                }
            } else // store hull points for handling later
                hull[d.triangles[i]] = c1;
        }
        // add hull connections in the form of "duplicated point" edges
        for (auto[i, c]: hull) {
            auto i1 = d.hull_tri[i];
            auto i2 = d.hull_tri[d.hull_next[i]];
            auto p1 = glm::vec2(d.coords[2 * d.triangles[i1]], d.coords[2 * d.triangles[i1] + 1]);
            auto p2 = glm::vec2(d.coords[2 * d.triangles[i2]], d.coords[2 * d.triangles[i2] + 1]);
            auto k1 = Diagram::site_key(p1);
            auto k2 = Diagram::site_key(p2);
            voroni.add_to_hull(added_faces[k1]);
            voroni.add_to_hull(added_faces[k2]);
            voroni.add_edge(c, c, added_faces[k1], added_faces[k2]);
        }

        return voroni;
    }

    Diagram DelaunatorAlgorithm::construct_clamped_voroni_diagram(const std::vector<double> &coords, double x, double y, double width, double height) {
        Diagram voroni;
        delaunator::Delaunator d(coords);
        std::unordered_map<std::string, unsigned int> added_faces;
        std::unordered_map<unsigned int, glm::vec2> hull;
        for (std::size_t i = 0; i < d.triangles.size(); i++) {
            auto j = d.halfedges[i];
            auto i0 = i - (i % 3);
            auto c1 = utils::math::compute_triangle_circumcenter(
                    glm::vec2(d.coords[2 * d.triangles[i0]], d.coords[2 * d.triangles[i0] + 1]),
                    glm::vec2(d.coords[2 * d.triangles[i0 + 1]], d.coords[2 * d.triangles[i0 + 1] + 1]),
                    glm::vec2(d.coords[2 * d.triangles[i0 + 2]], d.coords[2 * d.triangles[i0 + 2] + 1])
            );
            glm::vec2 site(d.coords[2 * d.triangles[i]], d.coords[2 * d.triangles[i] + 1]);
            auto k = Diagram::site_key(site);
            if (!added_faces.contains(k))
                added_faces[k] = voroni.add_face(site);
            if (j != delaunator::INVALID_INDEX) {
                if (i > j) {
                    auto j0 = j - (j % 3);
                    auto p2 = glm::vec2(d.coords[2 * d.triangles[j]], d.coords[2 * d.triangles[j] + 1]);
                    auto c2 = utils::math::compute_triangle_circumcenter(
                            glm::vec2(d.coords[2 * d.triangles[j0]], d.coords[2 * d.triangles[j0] + 1]),
                            glm::vec2(d.coords[2 * d.triangles[j0 + 1]], d.coords[2 * d.triangles[j0 + 1] + 1]),
                            glm::vec2(d.coords[2 * d.triangles[j0 + 2]], d.coords[2 * d.triangles[j0 + 2] + 1])
                    );
                    voroni.add_edge(c1, c2, added_faces[k], added_faces[Diagram::site_key(p2)]);
                }
            } else // store hull points for handling later
                hull[d.triangles[i]] = c1;
        }
        // add hull connections in the form of "duplicated point" edges
        for (auto[i, c]: hull) {
            auto i1 = d.hull_tri[i];
            auto i2 = d.hull_tri[d.hull_next[i]];
            auto p1 = glm::vec2(d.coords[2 * d.triangles[i1]], d.coords[2 * d.triangles[i1] + 1]);
            auto p2 = glm::vec2(d.coords[2 * d.triangles[i2]], d.coords[2 * d.triangles[i2] + 1]);
            auto k1 = Diagram::site_key(p1);
            auto k2 = Diagram::site_key(p2);
            voroni.add_to_hull(added_faces[k1]);
            voroni.add_to_hull(added_faces[k2]);
            voroni.add_edge(c, c, added_faces[k1], added_faces[k2]);
        }

        // cull edges completely outside of bounds first then modify edges that intersect with boundary
        auto max_x = x + width;
        auto max_y = y + height;
        for(auto &face: voroni.get_faces()) {
            auto neighbors = face.neighboring_edges;
            for (auto[n, e]: neighbors) {
                auto edge = voroni.get_edges()[e];
                // todo handle edges with endpoints out of bounds that still collide with boundary
                if (!collides(edge.first, x, y, max_x, max_y) && !collides(edge.second, x, y, max_x, max_y))
                    voroni.cull_edge(e);
            }
        }
        for(auto &face: voroni.get_faces()) {
            for(auto [n, e]: face.neighboring_edges) {
                auto edge = voroni.get_edges()[e];
                auto p1 = edge.first;
                auto p2 = edge.second;
                try {
                    if (!collides(p1, x, y, max_x, max_y)) {
                        if (p1.x < x) {
                            edge.second = compute_intersection(p1, p2, glm::vec2(x, y), glm::vec2(x, y + height));
                        } else if (p1.y < y) {
                            edge.second = compute_intersection(p1, p2, glm::vec2(x, y), glm::vec2(x + width, y));
                        } else if (p1.x >= max_x) {
                            edge.second = compute_intersection(p1, p2, glm::vec2(x + width, y),
                                                               glm::vec2(x + width, y + height));
                        } else if (p1.y >= max_y) {
                            edge.second = compute_intersection(p1, p2, glm::vec2(x, y + height),
                                                               glm::vec2(x + width, y + height));
                        }
                    } else if (!collides(p2, x, y, max_x, max_y)) {
                        if (p2.x < x) {
                            edge.first = compute_intersection(p1, p2, glm::vec2(x, y), glm::vec2(x, y + height));
                        } else if (p2.y < y) {
                            edge.first = compute_intersection(p1, p2, glm::vec2(x, y), glm::vec2(x + width, y));
                        } else if (p2.x >= max_x) {
                            edge.first = compute_intersection(p1, p2, glm::vec2(x + width, y),
                                                              glm::vec2(x + width, y + height));
                        } else if (p2.y >= max_y) {
                            edge.first = compute_intersection(p1, p2, glm::vec2(x, y + height),
                                                              glm::vec2(x + width, y + height));
                        }
                    }
                } catch(std::runtime_error ex) {
                    std::cout << e << ": " << ex.what() << std::endl;
                }
            }
        }

        return voroni;
    }

    Diagram DelaunatorAlgorithm::construct_delauney_diagram(const std::vector<double> &coords) {
        Diagram delauney;
        delaunator::Delaunator d(coords);
        std::unordered_map<std::string, unsigned int> added_faces;
        for (std::size_t i = 0; i < d.triangles.size(); i++) {
            auto i0 = i - (i % 3);
            auto site = utils::math::compute_triangle_circumcenter(
                    glm::vec2(d.coords[2 * d.triangles[i0]], d.coords[2 * d.triangles[i0] + 1]),
                    glm::vec2(d.coords[2 * d.triangles[i0 + 1]], d.coords[2 * d.triangles[i0 + 1] + 1]),
                    glm::vec2(d.coords[2 * d.triangles[i0 + 2]], d.coords[2 * d.triangles[i0 + 2] + 1])
            );
            if(!added_faces.contains(Diagram::site_key(site)))
                delauney.add_face(site);
        }
        for (std::size_t i = 0; i < d.triangles.size(); i++) {
            auto j = d.halfedges[i];
            if(j != delaunator::INVALID_INDEX) {
                glm::vec2 p1(d.coords[2 * d.triangles[i]], d.coords[2 * d.triangles[i] + 1]);
                glm::vec2 p2(d.coords[2 * d.triangles[j]], d.coords[2 * d.triangles[j] + 1]);
                delauney.add_edge(p1, p2, i, j);
            }
        }

        return delauney;
    }
} // namespace mapgen
