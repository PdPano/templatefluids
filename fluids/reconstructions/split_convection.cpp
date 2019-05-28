#include "split_convection.hpp"
#include "../derivatives/derivatives.hpp"
#include "../grid/cartesian_grid.hpp"
#include "../utils/flux_def.hpp"
#include "../utils/operators_overloads.hpp"
#include "../utils/point_def.hpp"
#include "../utils/point_functions.hpp"
#include "flux_functions/flux_interface.hpp"

#include <utility>

SplitConvection::SplitConvection(PointFunctions& pf_in,
    std::shared_ptr<Derivatives> der_in, std::shared_ptr<FluxInterface> flux_in)
    : Convection(pf_in, std::move(der_in))
    , flux(std::move(flux_in))
{
}

Flux SplitConvection::convection_x(const CartesianGrid& grid, int ind) const
{
    auto negative = [&](const Point& p) { return flux->fluxXNegative(p); };
    auto positive = [&](const Point& p) { return flux->fluxXPositive(p); };
    return -(der->DXForward(grid, negative, ind)
        + der->DXBackward(grid, positive, ind));
}

Flux SplitConvection::convection_y(const CartesianGrid& grid, int ind) const
{
    auto negative = [&](const Point& p) { return flux->fluxYNegative(p); };
    auto positive = [&](const Point& p) { return flux->fluxYPositive(p); };
    return -(der->DYForward(grid, negative, ind)
        + der->DYBackward(grid, positive, ind));
}
