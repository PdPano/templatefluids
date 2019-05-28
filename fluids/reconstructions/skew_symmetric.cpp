#include "skew_symmetric.hpp"
#include "../derivatives/derivatives.hpp"
#include "../grid/cartesian_grid.hpp"
#include "../utils/flux_def.hpp"
#include "../utils/point_def.hpp"
#include "../utils/point_functions.hpp"
#include "../utils/useful_alias.hpp"

#include <utility>

using alias::E;
using alias::P;
using alias::PointProperty;
using alias::RU;
using alias::RU2;
using alias::RUV;
using alias::RV;
using alias::RV2;
using alias::U;
using alias::V;

SkewSymmetric::SkewSymmetric(
    PointFunctions& pf_in, std::shared_ptr<Derivatives> der_in)
    : Convection(pf_in, std::move(der_in))
{
}

Flux SkewSymmetric::convection_x(const CartesianGrid& grid, int ind) const
{
    auto D = [&](PointProperty func) { return der->DX(grid, func, ind); };
    auto point = grid.values(ind);
    double ru = pf.ru(point);
    double u = pf.u(point);
    double v = pf.v(point);

    double d_rho = -D(RU);
    double d_ru = -1 / 2. * (D(RU2) + u * D(RU) + ru * D(U)) - D(P);
    double d_rv = -1 / 2. * (D(RUV) + v * D(RU) + ru * D(V));

    auto e_flux
        = [&](const Point& p) { return (pf.e(p) + pf.pressure(p)) * pf.u(p); };
    double full_e_conv = der->DX(grid, e_flux, ind);
    double d_e = -1 / 2. * (full_e_conv + (D(E) + D(P)) * u
                               + (pf.e(point) + pf.pressure(point)) * D(U));

    return {d_rho, d_ru, d_rv, d_e};
}

Flux SkewSymmetric::convection_y(const CartesianGrid& grid, int ind) const
{
    auto D = [&](PointProperty func) { return der->DY(grid, func, ind); };
    auto point = grid.values(ind);
    double rv = pf.rv(point);
    double u = pf.u(point);
    double v = pf.v(point);

    double d_rho = -D(RV);
    double d_ru = -1 / 2. * (D(RUV) + u * D(RV) + rv * D(U));
    double d_rv = -1 / 2. * (D(RV2) + v * D(RV) + rv * D(V)) - D(P);

    auto e_flux
        = [&](const Point& p) { return (pf.e(p) + pf.pressure(p)) * pf.v(p); };
    double full_e_conv = der->DY(grid, e_flux, ind);
    double d_e = -1 / 2. * (full_e_conv + (D(E) + D(P)) * v
                               + (pf.e(point) + pf.pressure(point)) * D(V));

    return {d_rho, d_ru, d_rv, d_e};
}
