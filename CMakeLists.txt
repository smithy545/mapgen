cmake_minimum_required(VERSION 3.16)

add_library(mapgen
        src/FortuneAlgorithm.cpp
        src/Terrain.cpp)
target_link_libraries(mapgen fmt utils engine glm)
target_include_directories(mapgen
        PRIVATE ${CMAKE_SOURCE_DIR}/external/delaunator-cpp/include
        PUBLIC include)
