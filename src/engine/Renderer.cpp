//
// Created by Philip Smith on 10/17/2020.
//

#include "Renderer.h"

#include <glm/gtx/rotate_vector.hpp>
#include <iostream>


bool Renderer::init(int width, int height) {
    assert(init_glfw());
    assert(init_window(width, height));
    assert(init_glew());
    assert(init_shaders());

    screen_width = width;
    screen_height = height;

    // cull triangles facing away from camera
    //glEnable(GL_CULL_FACE);

    // enable depth buffer
    //glEnable(GL_DEPTH_TEST);

    // background
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    return m_window;
}

bool Renderer::init_glfw() {
    if (!glfwInit())
        return false;

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    return true;
}

bool Renderer::init_window(int width, int height) {
    // Open a window and create its OpenGL context
    m_window = glfwCreateWindow(width, height, "Mapgen C++", nullptr, nullptr);
    if (m_window == nullptr) {
        std::cerr << "Failed to open GLFW window" << std::endl;
        glfwTerminate();
        return false;
    }
    glfwMakeContextCurrent(m_window);

    return true;
}

bool Renderer::init_glew() {
    glewExperimental = true; // Needed for core profile
    if (glewInit() != GLEW_OK) {
        glfwTerminate();
        return false;
    }

    return true;
}

bool Renderer::init_shaders() {
    const char *vshader_src = "#version 430 core\n"
                              "\n"
                              "layout (location = 0) in vec3 inPos;\n"
                              "layout (location = 1) in vec3 inColor;\n"
                              "layout (location = 2) in mat4 instanceModel;\n"
                              "\n"
                              "uniform mat4 VP;\n"
                              "\n"
                              "out vec4 color;\n"
                              "\n"
                              "void main() {\n"
                              "    gl_Position = VP * instanceModel * vec4(inPos, 1);\n"
                              "    color = vec4(inColor, 1);\n"
                              "}";
    int vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, &vshader_src, nullptr);
    glCompileShader(vertex_shader);
    int success;
    char infoLog[512];
    glGetShaderiv(vertex_shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(vertex_shader, 512, nullptr, infoLog);
        std::cout << "Error on vertex compilation: " << infoLog << std::endl;
        return false;
    }
    const char *fshader_src = "#version 430 core\n"
                              "\n"
                              "\n"
                              "in vec4 color;\n"
                              "out vec4 fragColor;\n"
                              "\n"
                              "\n"
                              "void main() {\n"
                              "    fragColor = color;\n"
                              "}";
    int fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, &fshader_src, nullptr);
    glCompileShader(fragment_shader);
    glGetShaderiv(fragment_shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(fragment_shader, 512, nullptr, infoLog);
        std::cout << "Error on fragment compilation: " << infoLog << std::endl;
        return false;
    }
    current_shader = glCreateProgram();
    glAttachShader(current_shader, vertex_shader);
    glAttachShader(current_shader, fragment_shader);
    glLinkProgram(current_shader);
    glGetProgramiv(current_shader, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(current_shader, 512, nullptr, infoLog);
        std::cout << "Error linking program: " << infoLog << std::endl;
        return false;
    }
    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
    glUseProgram(current_shader);

    return true;
}

void Renderer::load_scene(const Scene &scene) {
    std::size_t sizeof_vec4 = sizeof(glm::vec4);
    for (auto instance_list: scene.instances) {
        // setup vertex data
        GLuint vao;
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);

        // verts
        GLuint vbo;
        glGenBuffers(1, &vbo);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * instance_list.mesh.vertices.size(),
                     &instance_list.mesh.vertices[0], GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

        // colors
        GLuint cbo;
        glGenBuffers(1, &cbo);
        glBindBuffer(GL_ARRAY_BUFFER, cbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * instance_list.mesh.colors.size(),
                     &instance_list.mesh.colors[0], GL_STATIC_DRAW);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

        // indices
        GLuint ebo;
        glGenBuffers(1, &ebo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * instance_list.mesh.indices.size(),
                     &instance_list.mesh.indices[0], GL_STATIC_DRAW);

        // setup instance data
        GLuint ibo;
        glGenBuffers(1, &ibo);
        glBindBuffer(GL_ARRAY_BUFFER, ibo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(glm::mat4) * instance_list.instances.size(), &instance_list.instances[0],
                     GL_STATIC_DRAW);

        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 4 * sizeof_vec4, (void *) 0);
        glEnableVertexAttribArray(3);
        glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, 4 * sizeof_vec4, (void *) (1 * sizeof_vec4));
        glEnableVertexAttribArray(4);
        glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, 4 * sizeof_vec4, (void *) (2 * sizeof_vec4));
        glEnableVertexAttribArray(5);
        glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, 4 * sizeof_vec4, (void *) (3 * sizeof_vec4));

        glVertexAttribDivisor(2, 1);
        glVertexAttribDivisor(3, 1);
        glVertexAttribDivisor(4, 1);
        glVertexAttribDivisor(5, 1);

        vaos[instance_list.name] = {
                vao, // Vertex Array Object
                vbo, // Vertex Buffer Object
                cbo, // Color Buffer Object
                ebo, // Element Buffer Object
                instance_list.render_strategy, // GL_POINTS, GL_LINES or GL_TRIANGLES
                static_cast<unsigned int>(instance_list.mesh.indices.size()), // Number of indices
                static_cast<unsigned int>(instance_list.instances.size())};   // Number of instances
        glBindVertexArray(0);
    }
}

void Renderer::render_current_scene() {
    auto vp_uniform = glGetUniformLocation(current_shader, "VP");
    glm::mat4 view = glm::lookAt(camera.position,camera.position + camera.forward, camera.up);
    //glm::mat4 vp = glm::perspective(45.0f, screen_width / screen_height, 0.1f, 1000.0f) * view;
    glm::mat4 vp = glm::ortho(-0.f*screen_width, 1.f*screen_width, 1.f*screen_height, -0.f*screen_height);
    glUniformMatrix4fv(vp_uniform, 1, GL_FALSE, &vp[0][0]);
    for (const auto &vao: vaos) {
        glBindVertexArray(vao.second.vao);
        switch (vao.second.strategy) {
            case Scene::InstanceList::POINTS:
                glDrawElementsInstanced(GL_POINTS, vao.second.num_indices, GL_UNSIGNED_INT, 0, vao.second.num_instances);
                break;
            case Scene::InstanceList::LINES:
                glDrawElementsInstanced(GL_LINES, vao.second.num_indices, GL_UNSIGNED_INT, 0, vao.second.num_instances);
                break;
            case Scene::InstanceList::TRIANGLES:
                glDrawElementsInstanced(GL_TRIANGLES, vao.second.num_indices, GL_UNSIGNED_INT, 0, vao.second.num_instances);
                break;
            default:
                std::cerr << "Render strategy not supported" << std::endl;
        }
    }
}

void Renderer::resize(int width, int height) {
    screen_width = width;
    screen_height = height;
}

void Renderer::cleanup() {
    for (const auto &vao: vaos) {
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindVertexArray(0);

        glDeleteVertexArrays(1, &vao.second.vao);
        glDeleteBuffers(1, &vao.second.vbo);
        glDeleteBuffers(1, &vao.second.cbo);
        glDeleteBuffers(1, &vao.second.ebo);
    }
    glDeleteProgram(current_shader);
    glfwTerminate();
}

void Renderer::pan_horizontal(float diff) {
    camera.forward = glm::normalize(glm::rotate(camera.forward, 0.001f * diff, camera.up));
}

void Renderer::pan_vertical(float diff) {
    auto left = glm::cross(camera.up, camera.forward);
    camera.forward = glm::normalize(glm::rotate(camera.forward, 0.001f * diff, left));
}

void Renderer::move_forward() {
    camera.position.x += camera.scale * camera.forward.x;
    camera.position.z += camera.scale * camera.forward.z;
}

void Renderer::move_backward() {
    camera.position.x -= camera.scale * camera.forward.x;
    camera.position.z -= camera.scale * camera.forward.z;
}

void Renderer::move_left() {
    camera.position -= camera.scale * glm::cross(camera.forward, camera.up);
}

void Renderer::move_right() {
    camera.position += camera.scale * glm::cross(camera.forward, camera.up);
}

void Renderer::move_up() {
    camera.position += camera.scale * camera.up;
}

void Renderer::move_down() {
    camera.position -= camera.scale * camera.up;
}
