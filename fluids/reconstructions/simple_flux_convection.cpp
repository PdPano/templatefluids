#include "simple_flux_convection.hpp"
#include "../derivatives/derivatives.hpp"
#include "../grid/cartesian_grid.hpp"
#include "../utils/flux_def.hpp"
#include "../utils/operators_overloads.hpp"
#include "../utils/point_def.hpp"
#include "flux_functions/flux_interface.hpp"

#include <utility>

SimpleFluxConvection::SimpleFluxConvection(PointFunctions& pf_in,
    std::shared_ptr<Derivatives> der_in, std::shared_ptr<FluxInterface> flux_in)
    : Convection(pf_in, std::move(der_in))
    , flux(std::move(flux_in))
{
}

Flux SimpleFluxConvection::convection_x(
    const CartesianGrid& grid, int ind) const
{
    auto flux_func = [&](const Point& p) { return flux->fluxX(p); };
    return -der->DX(grid, flux_func, ind);
}

Flux SimpleFluxConvection::convection_y(
    const CartesianGrid& grid, int ind) const
{
    auto flux_func = [&](const Point& p) { return flux->fluxY(p); };
    return -der->DY(grid, flux_func, ind);
}
