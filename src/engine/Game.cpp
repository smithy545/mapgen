//
// Created by Philip Smith on 10/17/2020.
//

#include "Game.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <pathing/data.h>
#include <random>
#include <unordered_set>
#include <vector>

#include "pathing/FortuneAlgorithm.h"


static int screen_width;
static int screen_height;
static bool resized = false;
State Game::state;


void key_cb(GLFWwindow *window, int key, int scancode, int action, int mods) {
    switch (action) {
        case GLFW_PRESS:
            Game::state.set_key(key, true);
            break;
        case GLFW_RELEASE:
            Game::state.set_key(key, false);
            break;
        case GLFW_REPEAT:
            break;
        default:
            std::cerr << "Key action \"" << action << "\" not handled" << std::endl;
    }
}

void char_cb(GLFWwindow *window, unsigned int codepoint) {}

void cursor_pos_cb(GLFWwindow *window, double xpos, double ypos) {
    Game::state.set_mouse_x(xpos);
    Game::state.set_mouse_y(ypos);
    glfwSetCursorPos(window, screen_width / 2, screen_height / 2);
}

void mouse_scroll_cb(GLFWwindow *window, double xoffset, double yoffset) {
    Game::state.set_mouse_scroll(yoffset);
}

void mouse_button_cb(GLFWwindow *window, int button, int action, int mods) {}

void resize_cb(GLFWwindow *window, int width, int height) {
    glViewport(0, 0, width, height);
    screen_width = width;
    screen_height = height;
    resized = true;
}

void Game::run() {
    init();

    auto window = renderer.get_window();

    // setup cursor coords
    auto sw2 = screen_width / 2;
    auto sh2 = screen_height / 2;
    state.set_mouse_x(sw2);
    state.set_mouse_y(sh2);
    glfwSetCursorPos(window, sw2, sh2);

    state.start();
    state.pause();
    do {
        state.enter_frame();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        if (resized) {
            renderer.resize(screen_width, screen_height);
            resized = false;
        }

        if(!state.is_paused()) {
            if (!update())
                state.pause();
        }

        if (state.get_key(GLFW_KEY_W))
            renderer.move_forward();
        if (state.get_key(GLFW_KEY_A))
            renderer.move_left();
        if (state.get_key(GLFW_KEY_S))
            renderer.move_backward();
        if (state.get_key(GLFW_KEY_D))
            renderer.move_right();
        if(state.get_key(GLFW_KEY_SPACE))
            state.is_paused() ? state.unpause() : state.pause();
        if (state.get_key(GLFW_KEY_LEFT_SHIFT))
            renderer.move_up();
        if (state.get_key(GLFW_KEY_LEFT_CONTROL))
            renderer.move_down();
        if (state.get_mouse_x() != sw2)
            renderer.pan_horizontal(sw2 - state.get_mouse_x());
        if (state.get_mouse_y() != sh2)
            renderer.pan_vertical(state.get_mouse_y() - sh2);
        state.set_mouse_x(sw2);
        state.set_mouse_y(sh2);

        renderer.render_current_scene();

        glfwSwapBuffers(window);
        glfwPollEvents();
    } while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) == 0 &&
             !state.is_stopped());

    renderer.cleanup();
}

void Game::init() {
    screen_width = 1600;
    screen_height = 1200;
    assert(renderer.init(screen_width, screen_height));

    auto window = renderer.get_window();

    // Ensure we can capture the escape key being pressed below
    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    // cursor mode
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // Setup keyboard inputs
    glfwSetKeyCallback(window, key_cb);
    glfwSetCharCallback(window, char_cb);
    // Setup mouse inputs
    glfwSetCursorPosCallback(window, cursor_pos_cb);
    glfwSetScrollCallback(window, mouse_scroll_cb);
    glfwSetMouseButtonCallback(window, mouse_button_cb);
    // Window resize
    glfwSetFramebufferSizeCallback(window, resize_cb);

    init_world();
}

void Game::init_world() {
    std::vector<glm::mat4> instances{glm::mat4(1)};
    //std::vector<glm::vec2> sites;
    sites.clear();
    std::random_device rd;
    std::unordered_set<std::string> sitehashes;
    for (auto s: sites)
        sitehashes.insert(FortuneAlgorithm::site_key(s));

    for (int i = 0; i < 10000; i++) {
        glm::vec2 v;
        do
            v = glm::vec2(rd() % screen_width, rd() % screen_height);
        while (sitehashes.contains(FortuneAlgorithm::site_key(v)));
        sitehashes.insert(FortuneAlgorithm::site_key(v));
        //std::cout << FortuneAlgorithm::site_key(v) << ", ";
        sites.push_back(v);
    }
    std::sort(sites.begin(), sites.end(), FortuneAlgorithm::voroni_sort);
    std::cout << "\nInput sites generated." << std::endl;
    auto dt = std::chrono::system_clock::now().time_since_epoch();
    FortuneAlgorithm diagram(sites);
    auto voroni = diagram.construct(screen_width, screen_height, false);
    dt = std::chrono::system_clock::now().time_since_epoch() - dt;
    std::cout << "Generated voroni diagram from " << sites.size() << " sites in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(dt).count() << " ms with "
              << voroni.edges.size() << " edges" << std::endl;

    // site mesh
    std::vector<glm::vec3> site_verts;
    std::vector<glm::vec3> site_colors;
    std::vector<unsigned int> site_indices;
    // edge mesh
    std::vector<glm::vec3> verts;
    std::vector<glm::vec3> colors;
    std::vector<unsigned int> indices;
    for (auto site: sites) {
        site_verts.emplace_back(site.x, site.y, 0);
        site_colors.emplace_back(1, 0, 0);
        site_indices.push_back(site_indices.size());
    }
    for (auto edge: voroni.edges) {
        verts.emplace_back(edge.first.x, edge.first.y, 0);
        colors.emplace_back(0, 0, 1);
        indices.push_back(indices.size());
        verts.emplace_back(edge.second.x,  edge.second.y, 0);
        colors.emplace_back(0, 0, 1);
        indices.push_back(indices.size());
    }
    Scene scene;
    Scene::Mesh site_mesh{site_verts, site_colors, site_indices};
    Scene::Mesh edge_mesh{verts, colors, indices};
    scene.instances = {
            Scene::InstanceList{site_mesh, "sites", instances, Scene::InstanceList::POINTS},
            Scene::InstanceList{edge_mesh, "edges", instances, Scene::InstanceList::LINES}
    };
    std::cout << "Scene generated." << std::endl;
    renderer.load_scene(scene);
    std::cout << "Scene loaded." << std::endl;
}

bool Game::update() {
    std::random_device rd;
    std::vector<glm::mat4> instances{glm::mat4(1)};
    std::unordered_set<std::string> sitehashes;
    for (auto s: sites)
        sitehashes.insert(FortuneAlgorithm::site_key(s));
    for(auto & site : sites) {
        sitehashes.erase(FortuneAlgorithm::site_key(site));
        do {
            site.x = site.x + (rd() % 11) - 5;
            site.y = site.y + (rd() % 11) - 5;
        } while(sitehashes.contains(FortuneAlgorithm::site_key(site)));
        sitehashes.insert(FortuneAlgorithm::site_key(site));
    }
    FortuneAlgorithm diagram(sites);
    auto voroni = diagram.construct(false);
    // site mesh
    std::vector<glm::vec3> site_verts;
    std::vector<glm::vec3> site_colors;
    std::vector<unsigned int> site_indices;
    // edge mesh
    std::vector<glm::vec3> verts;
    std::vector<glm::vec3> colors;
    std::vector<unsigned int> indices;
    for (auto site: sites) {
        site_verts.emplace_back(site.x, site.y, 0);
        site_colors.emplace_back(1, 0, 0);
        site_indices.push_back(site_indices.size());
    }
    for (auto edge: voroni.edges) {
        verts.emplace_back(edge.first.x, edge.first.y, 0);
        colors.emplace_back(0, 0, 1);
        indices.push_back(indices.size());
        verts.emplace_back(edge.second.x, edge.second.y, 0);
        colors.emplace_back(0,0,1);
        indices.push_back(indices.size());
    }
    Scene scene;
    Scene::Mesh site_mesh{site_verts, site_colors, site_indices};
    Scene::Mesh edge_mesh{verts, colors, indices};
    scene.instances = {
            Scene::InstanceList{site_mesh, "sites", instances, Scene::InstanceList::POINTS},
            Scene::InstanceList{edge_mesh, "edges", instances, Scene::InstanceList::LINES}
    };
    renderer.load_scene(scene);
    return true;
}
