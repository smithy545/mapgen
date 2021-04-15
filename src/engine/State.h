//
// Created by Philip Smith on 10/18/2020.
//

#ifndef MAPGEN_STATE_H
#define MAPGEN_STATE_H


#include <chrono>
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "macros.h"


class State {
// game state
BVAR(paused, public, private){true};
BVAR(stopped, public, private){true};
VAR(double, fps, public, private){0};
VAR(std::chrono::time_point<std::chrono::system_clock>, last_frame_start, public, private);
// io state
VAR_GET(double, mouse_x, public){0};
VAR_GET(double, mouse_y, public){0};
VAR_GET(double, last_mouse_x, public){0};
VAR_GET(double, last_mouse_y, public){0};
VAR(double, mouse_scroll, public, public){0};
public:
    PTR(State);

    State();

    bool get_key(int key) const;

    void set_key(int key, bool value);

    void set_mouse_x(double x);

    void set_mouse_y(double y);

    void pause();

    void unpause();

    void start();

    void stop();

    double enter_frame();

private:
    bool keys[GLFW_KEY_LAST]{};
};


#endif //MAPGEN_STATE_H
