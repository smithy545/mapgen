cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 20)

add_library(mapgen
        src/terrain.cpp
        src/topology.cpp)
target_link_libraries(mapgen EnTT fmt utils engine BulletCollision BulletDynamics LinearMath)
target_link_directories(utils PRIVATE /opt/homebrew/Cellar/assimp/5.0.1/lib /opt/homebrew/Cellar/bullet/3.09/lib)
target_include_directories(mapgen
        PRIVATE
        ../../external/delaunator-cpp/include
        /opt/homebrew/Cellar/glm/0.9.9.8/include
        /opt/homebrew/Cellar/bullet/3.09/include/bullet
        PUBLIC include)
