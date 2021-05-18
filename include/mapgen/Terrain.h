//
// Created by Philip Smith on 4/17/2021.
//

#ifndef MAPGEN_TERRAIN_H
#define MAPGEN_TERRAIN_H

#include <entt/entt.hpp>
#include <mapgen/Diagram.h>
#include <string>
#include <unordered_set>
#include <utils/macros.h>


class Terrain {
public:
    Terrain(unsigned int num_sites, int width, int height, bool centered = false);

    entt::entity register_terrain_mesh(entt::registry& registry);

    void register_voroni_mesh(entt::registry& registry);
private:
    std::unordered_set<unsigned int> mountains;
    std::unordered_set<unsigned int> ocean;

    void assign_ocean(unsigned int start, int neighbor_depth);

    unsigned int find_nearest_mountain_face(unsigned int index);

    // TODO: add parameter to recursive search to stop searching after reaching neighbors 'X' hops away from index
    unsigned int find_nearest_mountain_face_recursive(unsigned int index);

    unsigned int find_mountain_kernel(unsigned int index, const std::unordered_set<unsigned int>& to_search, std::unordered_set<unsigned int> searched);

// underlying voroni/delauney diagrams
VAR_GET(Diagram, base, public);
VAR_GET(Diagram, dual, public);
};

#endif //MAPGEN_TERRAIN_H
