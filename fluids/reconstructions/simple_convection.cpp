#include "simple_convection.hpp"
#include "../derivatives/derivatives.hpp"
#include "../grid/cartesian_grid.hpp"
#include "../utils/flux_def.hpp"
#include "../utils/point_def.hpp"
#include "../utils/point_functions.hpp"

#include <utility>
using alias::P;
using alias::PointProperty;
using alias::RU;
using alias::RU2;
using alias::RUV;
using alias::RV;
using alias::RV2;

SimpleConvection::SimpleConvection(
    PointFunctions& pf_in, std::shared_ptr<Derivatives> der_in)
    : Convection(pf_in, std::move(der_in))
{
}

Flux SimpleConvection::convection_x(const CartesianGrid& grid, int ind) const
{
    auto D = [&](PointProperty func) { return der->DX(grid, func, ind); };

    double d_rho = -D(RU);
    double d_ru = -D(RU2) - D(P);
    double d_rv = -D(RUV);

    auto e_flux
        = [&](const Point& p) { return (pf.e(p) + pf.pressure(p)) * pf.u(p); };
    double d_e = -der->DX(grid, e_flux, ind);

    return {d_rho, d_ru, d_rv, d_e};
}

Flux SimpleConvection::convection_y(const CartesianGrid& grid, int ind) const
{
    auto D = [&](PointProperty func) { return der->DY(grid, func, ind); };

    double d_rho = -D(RV);
    double d_ru = -D(RUV);
    double d_rv = -D(RV2) - D(P);

    auto e_flux
        = [&](const Point& p) { return (pf.e(p) + pf.pressure(p)) * pf.v(p); };
    double d_e = -der->DY(grid, e_flux, ind);

    return {d_rho, d_ru, d_rv, d_e};
}
