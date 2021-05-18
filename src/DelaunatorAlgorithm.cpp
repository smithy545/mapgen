//
// Created by Philip Smith on 4/19/2021.
//

#include <mapgen/DelaunatorAlgorithm.h>
#include <utils/math_util.h>

namespace mapgen {
    Diagram DelaunatorAlgorithm::construct_voroni(const std::vector<double> &coords) {
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
} // namespace mapgen
