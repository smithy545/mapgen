cmake_minimum_required(VERSION 3.16)
project(mapgen C CXX)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenGL REQUIRED)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(external/glfw) # https://www.glfw.org/
add_subdirectory(external/glew) # http://glew.sourceforge.net/
add_subdirectory(external/fmt)  # https://github.com/fmtlib/fmt
add_subdirectory(external/json) # https://github.com/nlohmann/json
add_subdirectory(external/json-schema-validator) # https://github.com/pboettch/json-schema-validator
add_subdirectory(external/utils) # https://github.com/smithy545/utils

include_directories(
        src
        external/glm
        external/glew/include
)
link_libraries(
        glfw
        OpenGL::GL
        libglew_static
        nlohmann_json_schema_validator
        fmt
        utils
)

add_executable(mapgen
        src/main.cpp
        src/macros.h
        src/pathing/FortuneAlgorithm.cpp
        src/pathing/FortuneAlgorithm.h
        src/engine/Game.cpp
        src/engine/Game.h
        src/engine/Renderer.cpp
        src/engine/Renderer.h
        src/engine/Renderer.cpp
        src/engine/Renderer.h
        src/engine/Scene.h
        src/engine/State.cpp
        src/engine/State.h
        src/pathing/data.h
        src/pathing/Diagram.cpp
        src/pathing/Diagram.h)
