//
// Created by Philip Smith on 4/17/2021.
//

#ifndef CIVILWAR_TERRAIN_H
#define CIVILWAR_TERRAIN_H

#include <engine/Scene.h>
#include <mapgen/Diagram.h>
#include <unordered_map>
#include <utils/macros.h>


class Terrain {
public:
    Terrain(unsigned int num_sites, int width, int height, bool centered = false);

    Scene::Mesh get_terrain_mesh();

    // underlying voroni/delauney diagrams
    VAR_GET(Diagram, base, public);
    VAR_GET(Diagram, dual, public);
};

#endif //CIVILWAR_TERRAIN_H
