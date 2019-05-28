#include "split_convection_cached.hpp"
#include "../derivatives/derivatives.hpp"
#include "../grid/cartesian_grid.hpp"
#include "../utils/flux_def.hpp"
#include "../utils/operators_overloads.hpp"
#include "../utils/point_def.hpp"
#include "../utils/point_functions.hpp"
#include "flux_functions/flux_interface.hpp"

#include <utility>

SplitConvectionCached::SplitConvectionCached(PointFunctions& pf_in,
    std::shared_ptr<Derivatives> der_in, std::shared_ptr<FluxInterface> flux_in)
    : Convection(pf_in, std::move(der_in))
    , flux(std::move(flux_in))
{
}

void SplitConvectionCached::init(const CartesianGrid& grid)
{
    if (int(fluxXPos.size()) != grid.nPointsTotal) {
        fluxXPos.resize(grid.nPointsTotal);
    }
    if (int(fluxXNeg.size()) != grid.nPointsTotal) {
        fluxXNeg.resize(grid.nPointsTotal);
    }
    if (int(fluxYPos.size()) != grid.nPointsTotal) {
        fluxYPos.resize(grid.nPointsTotal);
    }
    if (int(fluxYNeg.size()) != grid.nPointsTotal) {
        fluxYNeg.resize(grid.nPointsTotal);
    }
#ifndef DEBUG
#pragma omp parallel for
#endif
    for (int ind = 0; ind < grid.nPointsTotal; ind++) {
        auto& p = grid.values(ind);
        fluxXPos[ind] = flux->fluxXPositive(p);
        fluxXNeg[ind] = flux->fluxXNegative(p);
        fluxYPos[ind] = flux->fluxYPositive(p);
        fluxYNeg[ind] = flux->fluxYNegative(p);
    }
}

Flux SplitConvectionCached::convection_x(
    const CartesianGrid& grid, int ind) const
{
    auto negative = [&](int ind) { return fluxXNeg[ind]; };
    auto positive = [&](int ind) { return fluxXPos[ind]; };
    return -(der->DXForward(grid, negative, ind)
        + der->DXBackward(grid, positive, ind));
}

Flux SplitConvectionCached::convection_y(
    const CartesianGrid& grid, int ind) const
{
    auto negative = [&](int ind) { return fluxXNeg[ind]; };
    auto positive = [&](int ind) { return fluxXPos[ind]; };
    return -(der->DYForward(grid, negative, ind)
        + der->DYBackward(grid, positive, ind));
}
