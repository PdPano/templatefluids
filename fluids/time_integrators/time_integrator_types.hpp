#ifndef TIME_INTEGRATOR_TYPES_HPP
#define TIME_INTEGRATOR_TYPES_HPP

#include "../utils/flux_def.hpp"
#include <vector>

struct CartesianVariation {
    CartesianVariation(int nPointsTotal)
        : grid_variation(std::vector<Flux>(nPointsTotal))
    {
    }
    std::vector<Flux> grid_variation;
};

#endif /* TIME_INTEGRATOR_TYPES_HPP */
