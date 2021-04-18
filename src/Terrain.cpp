//
// Created by Philip Smith on 4/17/2021.
//

#include <mapgen/Terrain.h>

#include <cassert>
#include <mapgen/FortuneAlgorithm.h>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>


Terrain::Terrain(unsigned int num_sites, int width, int height, bool centered) {
    assert(num_sites > 0);
    assert(width > 0);
    assert(height > 0);
    std::random_device rd;
    std::vector<double> coords;
    std::unordered_set<std::string> site_hashes;
    for (int i = 0; i < num_sites; i++) {
        glm::vec2 v;
        do {
            if(centered)
                v = glm::vec2(rd() % width - width / 2, rd() % height - height / 2);
            else
                v = glm::vec2(rd() % width, rd() % height);
        } while (site_hashes.contains(Diagram::site_key(v)));
        site_hashes.insert(Diagram::site_key(v));
        coords.push_back(v.x);
        coords.push_back(v.y);
    }
    m_base = FortuneAlgorithm::construct(coords);
    m_dual = m_base.dual();

    // assign cell types

    // assign elevations

    // generate terrain mesh
    for (int i = 0; i < m_base.get_faces().size(); i++) {
        auto face = m_base.get_faces()[i];
        auto color = glm::vec3(0,0,1);
        for(auto [neighbor, edge]: face.neighboring_edges) {
            if(neighbor < i)
                continue;
            auto e = m_base.get_edges()[edge];
            mesh.vertices.emplace_back(e.first.x, e.first.y, 0);
            mesh.colors.push_back(color);
            mesh.indices.push_back(mesh.indices.size());
            mesh.vertices.emplace_back(e.second.x, e.second.y, 0);
            mesh.colors.push_back(color);
            mesh.indices.push_back(mesh.indices.size());
        }
    }
    for (int i = 0; i < m_dual.get_faces().size(); i++) {
        auto face = m_dual.get_faces()[i];
        auto color = glm::vec3(1, 1, 0);
        for(auto [neighbor, edge]: face.neighboring_edges) {
            if(neighbor < i)
                continue;
            auto e = m_dual.get_edges()[edge];
            mesh.vertices.emplace_back(e.first.x, e.first.y, 0);
            mesh.colors.push_back(color);
            mesh.indices.push_back(mesh.indices.size());
            mesh.vertices.emplace_back(e.second.x,  e.second.y, 0);
            mesh.colors.push_back(color);
            mesh.indices.push_back(mesh.indices.size());
        }
    }
}

Scene::Mesh Terrain::get_terrain_mesh() const {
    return mesh;
}

Scene::Mesh Terrain::get_site_mesh() const {
    Scene::Mesh site_mesh;
    for (auto & face: m_base.get_faces()) {
        auto color = glm::vec3(1,0,0);
        site_mesh.vertices.emplace_back(face.site.x, face.site.y, 0);
        site_mesh.colors.push_back(color);
        site_mesh.indices.push_back(site_mesh.indices.size());
    }
    return site_mesh;
}
