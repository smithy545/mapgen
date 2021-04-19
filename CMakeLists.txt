cmake_minimum_required(VERSION 3.16)

add_subdirectory(external/delaunator-cpp)
add_library(mapgen
        src/FortuneAlgorithm.cpp
        src/Terrain.cpp
        src/DelaunatorAlgorithm.cpp)
target_link_libraries(mapgen EnTT fmt utils engine glm)
target_include_directories(mapgen
        PRIVATE external/delaunator-cpp/include
        PUBLIC include)
