//
// Created by Philip Smith on 4/17/2021.
//

#include <mapgen/Terrain.h>

#include <cassert>
#include <limits>
#include <map>
#include <mapgen/DelaunatorAlgorithm.h>
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
    m_base = DelaunatorAlgorithm::construct_voroni(coords);
    m_dual = m_base.dual();

    // assign ocean tiles
    for(auto index: m_base.get_hull())
        assign_ocean(index, 6);
}

entt::entity Terrain::register_terrain_mesh(entt::registry& registry, std::string id) {
    // setup hash table of face indices based on site hash for convenient access
    std::unordered_map<std::string, unsigned int> facii;
    for(int i = 0; i < m_base.get_faces().size(); i++) {
        auto face = m_base.get_faces()[i];
        facii[Diagram::site_key(face.site)] = i;
    }
    // assign elevations based on distance from ocean
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
        level *= 1.2;
    }
    // normalize elevation and assign mountain peaks based on elevation "probability"
    for(auto [i, v]: elevations) {
        elevations[i] = v/level;
        auto r = ((double) rand()/(RAND_MAX));
        if(r < elevations[i])
            mountains.insert(i);
    }

    for (int index = 0; index < m_base.get_faces().size(); index++) {
        if(!ocean.contains(index) && !mountains.contains(index)) {
            auto face = m_base.get_faces()[index];
            auto nearest = recursively_find_nearest_mountain_face(index);
            if (nearest != index) {
                auto d = glm::distance(m_base.get_faces()[nearest].site, face.site);
                // scale height inverse logarithmically to distance from mountains (range 0-1 when shifting x by e)
                if(d > 5.0)
                    elevations[index] *= 2.0/glm::log(d + glm::exp(1));
                else
                    elevations[index] = elevations[nearest]; // turn into plateau if close enough
            }
        }
    }

    // generate meshes
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> colors;
    std::vector<unsigned int> indices;
    std::unordered_set<std::string, unsigned int> added_points;
    for (int index = 0; index < m_base.get_faces().size(); index++) {
        auto face = m_base.get_faces()[index];
        float z = elevations[index];
        if(!ocean.contains(index))
            z = 400*elevations[index];
        auto ocean_color = glm::vec3(0, 0, 1);
        auto land_color  = glm::vec3(.5, .3, .1);
        auto mountain_color = glm::vec3(.6,.6,.6);
        auto mountain_top_color = glm::vec3(1, 1, 1);
        auto color = ocean_color;
        if(z > 0)
            color = land_color + glm::vec3(.2, .4, .4)/z;
        if(z > 100)
            color = mountain_color;
        if(z > 200)
            color = mountain_top_color;
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

unsigned int Terrain::find_nearest_mountain_face(unsigned int index) {
    auto distance = std::numeric_limits<double>::infinity();
    auto face = m_base.get_faces()[index];
    auto found = index;
    for(auto i: mountains) {
        if(index != i) {
            auto mountain_face = m_base.get_faces()[i];
            auto d = glm::distance(face.site, mountain_face.site);
            if(found == index || d < distance) {
                found = i;
                distance = d;
            }
        }
    }
    return found;
}

unsigned int Terrain::recursively_find_nearest_mountain_face(unsigned int index) {
    auto face = m_base.get_faces()[index];
    auto neighbors = face.neighboring_edges;
    std::unordered_set<unsigned int> next{};
    std::unordered_set<unsigned int> searched{index};
    auto distance = std::numeric_limits<double>::infinity();
    auto found = index;
    for (auto [n, e]: neighbors) {
        if (mountains.contains(n) && glm::distance(face.site, m_base.get_faces()[n].site) < distance) {
            found = n;
            distance = glm::distance(face.site, m_base.get_faces()[n].site);
        }
        searched.insert(n);
        for(auto [n2, e2]: m_base.get_faces()[n].neighboring_edges){
            if(!searched.contains(n2))
                next.insert(n2);
        }
    }
    if(found == index)
        return find_mountain_kernel(index, next, searched);
    return found;
}

unsigned int Terrain::find_mountain_kernel(unsigned int index, const std::unordered_set<unsigned int>& to_search, std::unordered_set<unsigned int> searched) {
    if(searched.size() == m_base.get_faces().size())
        return index;
    std::unordered_set<unsigned int> next{};
    auto distance = std::numeric_limits<double>::infinity();
    auto found = index;
    auto base = m_base.get_faces()[index];
    for (auto i: to_search) {
        auto face = m_base.get_faces()[i];
        if (mountains.contains(i) && glm::distance(face.site, base.site) < distance) {
            found = i;
            distance = glm::distance(face.site, base.site);
        }
        searched.insert(i);
        for(auto [n, e]: face.neighboring_edges){
            if(!searched.contains(n))
                next.insert(n);
        }
    }
    if(found == index)
        find_mountain_kernel(index, next, searched);
    return found;
}
