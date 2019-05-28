#include "simple_dissipation.hpp"

#include "../derivatives/derivatives.hpp"
#include "../grid/cartesian_grid.hpp"
#include "../utils/flux_def.hpp"
#include "../utils/point_def.hpp"
#include "../utils/point_functions.hpp"
#include <utility>
using alias::PointProperty;
using alias::U;
using alias::V;

SimpleDissipation::SimpleDissipation(PointFunctions& pf_in,
    std::shared_ptr<Derivatives> der_in, double reynolds_in, double prandtl_in)
    : Dissipation(pf_in, std::move(der_in), reynolds_in, prandtl_in)
{
}

Flux SimpleDissipation::dissipation_x(const CartesianGrid& grid, int ind) const
{

    double Txx = dissipation_tool.Txx(grid, ind);
    double Txy = dissipation_tool.Txy(grid, ind);
    double dTxx_dx = dissipation_tool.dTxx_dx(grid, ind);
    double dTxy_dx = dissipation_tool.dTxy_dx(grid, ind);
    double dqx_dx = dissipation_tool.dqx_dx(grid, ind);

    auto point = grid.values(ind);
    double e_diss = der->DX(grid, alias::U, ind) * Txx + pf.u(point) * dTxx_dx
        + der->DX(grid, alias::V, ind) * Txy + pf.v(point) * dTxy_dx + dqx_dx;
    return {0., dTxx_dx, dTxy_dx, e_diss};
}

Flux SimpleDissipation::dissipation_y(const CartesianGrid& grid, int ind) const
{
    double Txy = dissipation_tool.Txy(grid, ind);
    double Tyy = dissipation_tool.Tyy(grid, ind);
    double dTxy_dy = dissipation_tool.dTxy_dy(grid, ind);
    double dTyy_dy = dissipation_tool.dTyy_dy(grid, ind);
    double dqy_dy = dissipation_tool.dqy_dy(grid, ind);

    auto point = grid.values(ind);
    double e_diss = der->DY(grid, alias::V, ind) * Tyy + pf.v(point) * dTyy_dy
        + der->DY(grid, alias::U, ind) * Txy + pf.u(point) * dTxy_dy + dqy_dy;
    return {0., dTxy_dy, dTyy_dy, e_diss};
}
