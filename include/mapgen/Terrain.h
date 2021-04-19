//
// Created by Philip Smith on 4/17/2021.
//

#ifndef CIVILWAR_TERRAIN_H
#define CIVILWAR_TERRAIN_H

#include <engine/Scene.h>
#include <mapgen/Diagram.h>
#include <unordered_map>
#include <unordered_set>
#include <utils/macros.h>


class Terrain {
public:
    Terrain(unsigned int num_sites, int width, int height, bool centered = false);

    Scene::Mesh get_terrain_mesh() const;

    Scene::Mesh get_site_mesh() const;

    Scene::Mesh get_wireframe_mesh() const;
private:
    Scene::Mesh wire_mesh;
    Scene::Mesh terrain_mesh;
    std::unordered_set<unsigned int> ocean;
    std::unordered_set<unsigned int> beach;
    std::unordered_set<unsigned int> land;
    std::unordered_map<unsigned int, double> elevations;

    void assign_ocean(unsigned int start, int neighbor_depth);

// underlying voroni/delauney diagrams
VAR_GET(Diagram, base, public);
VAR_GET(Diagram, dual, public);
};

#endif //CIVILWAR_TERRAIN_H
