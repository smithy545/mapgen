//
// Created by Philip Smith on 10/17/2020.
//

#ifndef MAPGEN_SCENE_H
#define MAPGEN_SCENE_H

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <gl/GL.h>
#include <string>
#include <vector>


struct Scene {
    struct Mesh {
        std::vector<glm::vec3> vertices;
        std::vector<glm::vec3> colors;
        std::vector<unsigned int> indices;
    };
    struct InstanceList {
        Mesh mesh;
        std::string name;
        std::vector<glm::mat4> instances;
        enum RenderStrategy {
            POINTS,
            LINES,
            TRIANGLES
        } render_strategy{TRIANGLES};
    };
    struct Light {
        glm::vec3 position;
        double radius;
    };
    struct DirectedLight : Light {
        glm::vec3 direction;
    };
    std::vector<Light *> lights;
    std::vector<InstanceList> instances;
};


#endif //MAPGEN_SCENE_H
