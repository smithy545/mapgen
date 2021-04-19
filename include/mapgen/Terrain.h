//
// Created by Philip Smith on 4/17/2021.
//

#ifndef CIVILWAR_TERRAIN_H
#define CIVILWAR_TERRAIN_H

#include <engine/Scene.h>
#include <entt/entt.hpp>
#include <mapgen/Diagram.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utils/macros.h>


class Terrain {
public:
    Terrain(unsigned int num_sites, int width, int height, bool centered = false);

    entt::entity register_terrain_mesh(entt::registry& registry, std::string id);

    entt::entity register_site_mesh(entt::registry& registry, std::string id);

    entt::entity register_wireframe_mesh(entt::registry& registry, std::string id);
private:
    std::unordered_set<unsigned int> ocean;

    void assign_ocean(unsigned int start, int neighbor_depth);

// underlying voroni/delauney diagrams
VAR_GET(Diagram, base, public);
VAR_GET(Diagram, dual, public);
};

#endif //CIVILWAR_TERRAIN_H
