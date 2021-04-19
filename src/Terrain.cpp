//
// Created by Philip Smith on 4/17/2021.
//

#include <mapgen/Terrain.h>

#include <cassert>
#include <map>
#include <mapgen/FortuneAlgorithm.h>
#include <random>
#include <string>
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

    // assign ocean tiles
    for(auto index: m_base.get_hull())
        assign_ocean(index, 3);

    std::unordered_map<std::string, unsigned int> facii;
    for(unsigned int i = 0; i < m_base.get_faces().size(); i++) {
        auto face = m_base.get_faces()[i];
        facii[Diagram::site_key(face.site)] = i;
    }
    // assign elevations and tile type based on distance from ocean
    std::unordered_set<unsigned int> layer;
    for(auto index: ocean) {
        elevations[index] = 0;
        for(auto [n, e]: m_base.get_faces()[index].neighboring_edges) {
            if(!ocean.contains(n)) {
                layer.insert(n);
                beach.insert(n);
            }
        }
    }
    double level = 1.0;
    while(!layer.empty()) {
        std::unordered_set<unsigned int> next_layer;
        for(auto i: layer) {
            elevations[i] = level;
            for(auto [n, e]: m_base.get_faces()[i].neighboring_edges) {
                if(!elevations.contains(n)) {
                    next_layer.insert(n);
                    land.insert(n);
                }
            }
        }
        layer = next_layer;
        level += 10.0;
    }

    // generate meshes
    for (auto & edge: m_base.get_edges()) {
        auto color = glm::vec3(1, 1, 1);
        wire_mesh.vertices.emplace_back(edge.first.x, edge.first.y, 0);
        wire_mesh.vertices.emplace_back(edge.second.x, edge.second.y, 0);
        wire_mesh.colors.push_back(color);
        wire_mesh.colors.push_back(color);
        wire_mesh.indices.push_back(wire_mesh.indices.size());
        wire_mesh.indices.push_back(wire_mesh.indices.size());
    }
    for (auto & edge: m_dual.get_edges()) {
        auto color = glm::vec3(1, 1, 0);
        wire_mesh.vertices.emplace_back(edge.first.x, edge.first.y, 0);
        wire_mesh.vertices.emplace_back(edge.second.x, edge.second.y, 0);
        wire_mesh.colors.push_back(color);
        wire_mesh.colors.push_back(color);
        wire_mesh.indices.push_back(wire_mesh.indices.size());
        wire_mesh.indices.push_back(wire_mesh.indices.size());
    }
    std::unordered_set<std::string, unsigned int> added_points;
    for (int index = 0; index < m_base.get_faces().size(); index++) {
        auto face = m_base.get_faces()[index];

        auto ocean_color = glm::vec3(.0, .0, 1.);
        auto beach_color = glm::vec3(.7, .6, .3);
        auto land_color  = glm::vec3(.6, .3, .1);
        auto color = ocean_color;
        if(beach.contains(index))
            color = beach_color;
        else if(land.contains(index))
            color = land_color;
        terrain_mesh.vertices.emplace_back(face.site.x, elevations[index], face.site.y);
        terrain_mesh.colors.push_back(color);

        std::map<double, unsigned int> neighbors;
        for(auto [n, e]: face.neighboring_edges) {
            auto neighbor = m_base.get_faces()[n];
            neighbors[std::atan2(neighbor.site.y - face.site.y, neighbor.site.x - face.site.x)] = n;
        }
        auto last_angle = glm::two_pi<double>();
        for(auto [angle, n1]: neighbors) {
            if(last_angle != glm::two_pi<double>()) {
                auto n2 = neighbors[last_angle];
                terrain_mesh.indices.push_back(index);
                terrain_mesh.indices.push_back(n1);
                terrain_mesh.indices.push_back(n2);
            }
            last_angle = angle;
        }
    }
}

Scene::Mesh Terrain::get_terrain_mesh() const {
    return terrain_mesh;
}

Scene::Mesh Terrain::get_site_mesh() const {
    Scene::Mesh site_mesh;
    for (auto & face: m_base.get_faces()) {
        auto color = glm::vec3(1,0,0);
        site_mesh.vertices.emplace_back(face.site.x, 0, face.site.y);
        site_mesh.colors.push_back(color);
        site_mesh.indices.push_back(site_mesh.indices.size());
    }
    return site_mesh;
}

Scene::Mesh Terrain::get_wireframe_mesh() const {
    return wire_mesh;
}

void Terrain::assign_ocean(unsigned int start, int neighbor_depth) {
    if(neighbor_depth == 0)
        return;
    ocean.insert(start);
    for(auto [neighbor, edge]: m_base.get_faces()[start].neighboring_edges)
        assign_ocean(neighbor, neighbor_depth - 1);
}
