//
// Created by Philip Smith on 4/17/2021.
//

#include <mapgen/Terrain.h>

#include <cassert>
#include <mapgen/FortuneAlgorithm.h>
#include <random>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>


Terrain::Terrain(unsigned int num_sites, int width, int height, bool centered) {
    assert(num_sites > 0);
    assert(width > 0);
    assert(height > 0);
    std::random_device rd;
    std::vector<glm::vec2> sites;
    std::unordered_set<std::string> sitehashes;
    for (int i = 0; i < num_sites; i++) {
        glm::vec2 v;
        do {
            if(centered)
                v = glm::vec2(rd() % width - width / 2, rd() % height - height / 2);
            else
                v = glm::vec2(rd() % width, rd() % height);
        } while (sitehashes.contains(Diagram::site_key(v)));
        sitehashes.insert(Diagram::site_key(v));
        sites.push_back(v);
    }
    m_base = FortuneAlgorithm::construct_via_delaunator(sites);
    m_dual = m_base.dual();

    // assign elevations

}

Scene::Mesh Terrain::get_terrain_mesh() {
    Scene::Mesh terrain_mesh;

    // edge mesh
    for (auto & edge: m_base.get_edges()) {
        terrain_mesh.vertices.emplace_back(edge.first.x, edge.first.y, 0);
        terrain_mesh.colors.emplace_back(0, 0, 1);
        terrain_mesh.indices.push_back(terrain_mesh.indices.size());
        terrain_mesh.vertices.emplace_back(edge.second.x,  edge.second.y, 0);
        terrain_mesh.colors.emplace_back(0, 0, 1);
        terrain_mesh.indices.push_back(terrain_mesh.indices.size());
    }
    for (auto & edge: m_dual.get_edges()) {
        terrain_mesh.vertices.emplace_back(edge.first.x, edge.first.y, 0);
        terrain_mesh.colors.emplace_back(1, 0, 1);
        terrain_mesh.indices.push_back(terrain_mesh.indices.size());
        terrain_mesh.vertices.emplace_back(edge.second.x,  edge.second.y, 0);
        terrain_mesh.colors.emplace_back(1, 0, 1);
        terrain_mesh.indices.push_back(terrain_mesh.indices.size());
    }

    return terrain_mesh;
}
