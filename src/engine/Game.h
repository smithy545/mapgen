//
// Created by Philip Smith on 10/17/2020.
//

#ifndef MAPGEN_GAME_H
#define MAPGEN_GAME_H

#include "Renderer.h"
#include "State.h"


class Game {
public:
    static State state;

    void operator()() {
        run();
    }

    void run();

private:
    Renderer renderer{};

    void init();
    void init_world();
    bool update();
};


#endif //MAPGEN_GAME_H
