#include <iostream>
#include <thread>

#include "engine/Game.h"


int main() {
    Game game_obj{};
    std::thread game_thread(game_obj);
    game_thread.join();

    return 0;
}
