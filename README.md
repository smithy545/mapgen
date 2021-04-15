# Mapgen C++

WIP

External libs not included but references to source websites are in CMakeLists.txt. Currently built in Windows.

A C++ version of the redblobgames Mapgen4 (https://www.redblobgames.com/maps/mapgen4/)

Fortune's algorithm implementation heavily on the blogs and code samples published by Pierre Vigier:

https://pvigier.github.io/2018/11/18/fortune-algorithm-details.html

https://github.com/pvigier/FortuneAlgorithm

and Jaque Huenis:

https://jacquesheunis.com/post/fortunes-algorithm/

https://github.com/jacquesh/fortunes-algorithm

It extends on these implementations with "parallel arc line" handling (three adjacent arcs with equal y-coords) and "seam"
handling (an event occuring on the intersection of two existing arcs, essentially a combined circle/site event).

To-do:
- Bound voroni diagram
- Optimize memory usage by reducing dynamic memory allocation for static indexed memory
- Transition to smart pointers
- Write function to determine water cells
- Map diagram to elevations
- Write river pathing
- Write vegetation generation
- Wind?
- Code cleanup (lol)