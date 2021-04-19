//
// Created by Philip Smith on 4/17/2021.
//

#include <mapgen/Terrain.h>

#include <cassert>
#include <map>
#include <mapgen/FortuneAlgorithm.h>
#include <random>
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
}

entt::entity Terrain::register_terrain_mesh(entt::registry& registry, std::string id) {
    std::unordered_map<std::string, unsigned int> facii;
    for(int i = 0; i < m_base.get_faces().size(); i++) {
        auto face = m_base.get_faces()[i];
        facii[Diagram::site_key(face.site)] = i;
    }
    // assign elevations and tile type based on distance from ocean
    std::unordered_set<unsigned int> layer;
    std::unordered_map<unsigned int, double> elevations;
    for(auto index: ocean) {
        elevations[index] = 0;
        for(auto [n, e]: m_base.get_faces()[index].neighboring_edges) {
            if(!ocean.contains(n))
                layer.insert(n);
        }
    }
    double level = 1.0;
    while(!layer.empty()) {
        std::unordered_set<unsigned int> next_layer;
        for(auto i: layer) {
            elevations[i] = level;
            for(auto [n, e]: m_base.get_faces()[i].neighboring_edges) {
                if(!elevations.contains(n))
                    next_layer.insert(n);
            }
        }
        layer = next_layer;
        level += 10.0;
    }
    // normalize elevation
    auto A = 250.0;
    for(auto [i, v]: elevations) {
        elevations[i] = A*(v/level);
    }

    // generate meshes
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> colors;
    std::vector<unsigned int> indices;
    std::unordered_set<std::string, unsigned int> added_points;
    for (int index = 0; index < m_base.get_faces().size(); index++) {
        auto face = m_base.get_faces()[index];
        auto z = elevations[index];
        auto ocean_color = glm::vec3(0, 0, 1);
        auto beach_color = glm::vec3(.7, .6, .3);
        auto land_color  = glm::vec3(.6, .3, .1);
        auto color = ocean_color;
        if(z > 0)
            color = land_color;
        if(z > 100)
            color = glm::vec3(.6,.6,.6);
        if(z > 200)
            color = glm::vec3(1,1,1);
        vertices.emplace_back(face.site.x, z, face.site.y);
        colors.push_back(color);

        std::map<double, unsigned int> neighbors;
        for(auto [n, e]: face.neighboring_edges) {
            auto neighbor = m_base.get_faces()[n];
            neighbors[std::atan2(neighbor.site.y - face.site.y, neighbor.site.x - face.site.x)] = n;
        }
        auto last_angle = glm::two_pi<double>();
        for(auto [angle, n1]: neighbors) {
            if(last_angle != glm::two_pi<double>()) {
                auto n2 = neighbors[last_angle];
                indices.push_back(index);
                indices.push_back(n1);
                indices.push_back(n2);
            }
            last_angle = angle;
        }
    }
    auto mesh_entity = registry.create();
    registry.emplace<Scene::Mesh>(mesh_entity, vertices, colors, indices);
    registry.emplace<std::string>(mesh_entity, id);
    registry.emplace<Scene::InstanceList::RenderStrategy>(mesh_entity, Scene::InstanceList::TRIANGLES);
    return mesh_entity;
}

entt::entity Terrain::register_site_mesh(entt::registry& registry, std::string id) {
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> colors;
    std::vector<unsigned int> indices;
    for (auto & face: m_base.get_faces()) {
        auto color = glm::vec3(1,0,0);
        vertices.emplace_back(face.site.x, face.site.y, 0);
        colors.push_back(color);
        indices.push_back(indices.size());
    }
    auto mesh_entity = registry.create();
    registry.emplace<Scene::Mesh>(mesh_entity, vertices, colors, indices);
    registry.emplace<std::string>(mesh_entity, id);
    registry.emplace<Scene::InstanceList::RenderStrategy>(mesh_entity, Scene::InstanceList::POINTS);
    return mesh_entity;
}

entt::entity Terrain::register_wireframe_mesh(entt::registry& registry, std::string id) {
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> colors;
    std::vector<unsigned int> indices;
    for (auto & edge: m_base.get_edges()) {
        auto color = glm::vec3(1, 1, 1);
        vertices.emplace_back(edge.first.x, edge.first.y, 0);
        vertices.emplace_back(edge.second.x, edge.second.y, 0);
        colors.push_back(color);
        colors.push_back(color);
        indices.push_back(indices.size());
        indices.push_back(indices.size());
    }
    for (auto & edge: m_dual.get_edges()) {
        auto color = glm::vec3(1, 1, 0);
        vertices.emplace_back(edge.first.x, edge.first.y, 0);
        vertices.emplace_back(edge.second.x, edge.second.y, 0);
        colors.push_back(color);
        colors.push_back(color);
        indices.push_back(indices.size());
        indices.push_back(indices.size());
    }
    auto mesh_entity = registry.create();
    registry.emplace<Scene::Mesh>(mesh_entity, vertices, colors, indices);
    registry.emplace<std::string>(mesh_entity, id);
    registry.emplace<Scene::InstanceList::RenderStrategy>(mesh_entity, Scene::InstanceList::LINES);
    return mesh_entity;
}

void Terrain::assign_ocean(unsigned int start, int neighbor_depth) {
    if(neighbor_depth == 0)
        return;
    ocean.insert(start);
    for(auto [neighbor, edge]: m_base.get_faces()[start].neighboring_edges)
        assign_ocean(neighbor, neighbor_depth - 1);
}
