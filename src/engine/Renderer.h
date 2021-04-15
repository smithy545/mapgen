//
// Created by Philip Smith on 10/17/2020.
//

#ifndef MAPGEN_RENDERER_H
#define MAPGEN_RENDERER_H

#include <glm/ext.hpp>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <string>
#include <unordered_map>

#include "macros.h"
#include "Scene.h"


class Renderer {
VAR_GET(GLFWwindow*, window, public);
public:
    struct Camera {
        glm::vec3 position{0, 0, 0};
        glm::vec3 forward{0, 0, -1};
        glm::vec3 up{0, 1, 0};
        float scale{1.f};
    };

    struct VertexArrayObject {
        GLuint vao;
        GLuint vbo;
        GLuint cbo;
        GLuint ebo;
        Scene::InstanceList::RenderStrategy strategy;
        unsigned int num_indices;
        unsigned int num_instances;
    };

    bool init(int width, int height);

    void render_current_scene();

    void load_scene(const Scene &scene);

    void resize(int width, int height);

    void cleanup();

    void pan_horizontal(float diff);

    void pan_vertical(float diff);

    void move_forward();

    void move_backward();

    void move_left();

    void move_right();

    void move_up();

    void move_down();

private:
    float screen_width, screen_height;
    Camera camera;
    GLuint current_shader;
    std::unordered_map<std::string, VertexArrayObject> vaos;

    bool init_window(int width, int height);

    bool init_glfw();

    bool init_glew();

    bool init_shaders();
};


#endif //MAPGEN_RENDERER_H
