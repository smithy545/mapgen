//
// Created by Philip Smith on 4/19/2021.
//

#ifndef MAPGEN_DELAUNATORALGORITHM_H
#define MAPGEN_DELAUNATORALGORITHM_H

#include <delaunator-header-only.hpp>
#include <vector>

#include "Diagram.h"

namespace mapgen {
    class DelaunatorAlgorithm {
    public:
        static Diagram construct_voroni(const std::vector<double> &coords);
    };
} // namespace mapgen

#endif //MAPGEN_DELAUNATORALGORITHM_H
