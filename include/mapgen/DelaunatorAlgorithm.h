//
// Created by Philip Smith on 4/19/2021.
//

#ifndef CIVILWAR_DELAUNATORALGORITHM_H
#define CIVILWAR_DELAUNATORALGORITHM_H

#include <delaunator-header-only.hpp>
#include <vector>

#include "Diagram.h"


class DelaunatorAlgorithm {
public:
    static Diagram construct_voroni(const std::vector<double>& coords);
};

#endif //CIVILWAR_DELAUNATORALGORITHM_H
