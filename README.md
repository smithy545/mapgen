# Mapgen C++

WIP

A C++ version of the redblobgames Mapgen4 (https://www.redblobgames.com/maps/mapgen4/)

#### Installation
Libs used:

https://www.glfw.org/

http://glew.sourceforge.net/

https://github.com/fmtlib/fmt

https://github.com/nlohmann/json

https://github.com/pboettch/json-schema-validator

https://github.com/smithy545/utils

https://github.com/abellgithub/delaunator-cpp

#### Todo
- Transition to generating "mountain peaks" over generic elevation
- Write river pathing
- Write vegetation generation
- Wind?
- Code cleanup (lol)

#### Done so far

- Can generate the dual of a given diagram (Voroni Diagram -> Delauney Triangulation -> Voroni etc.)

- Generates basic terrain mesh 

- Basic renderer/input handling/game loop for displaying and updating Voroni Diagram in 2d (ortho) or 3d (perspective)

    - Renderer built on OpenGl
    - OpenGlMathematics (GLM) used for most math functions
    - Window and input system built on GLFW
    - Fmt string lib for string formatting
    - nlohmann and pboettch json libraries for json handling
    - Basic mesh system outline with EnTT for memory/event handling
