cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 20)

add_subdirectory(external/delaunator-cpp)
add_library(mapgen
        src/FortuneAlgorithm.cpp
        src/Terrain.cpp
        src/DelaunatorAlgorithm.cpp)
target_link_libraries(mapgen EnTT fmt utils engine)
target_include_directories(mapgen
        PRIVATE
        external/delaunator-cpp/include
        /opt/homebrew/Cellar/glm/0.9.9.8/include
        PUBLIC include)
