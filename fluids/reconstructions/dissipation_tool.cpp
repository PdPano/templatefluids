#include "dissipation_tool.hpp"
#include "../derivatives/derivatives.hpp"
#include "../grid/cartesian_grid.hpp"
#include "../utils/point_functions.hpp"
#include "../utils/useful_alias.hpp"

#include <utility>

DissipationTool::DissipationTool(const PointFunctions& pf_in,
    std::shared_ptr<Derivatives> der_in, double reynolds_in, double prandtl_in)
    : der(std::move(der_in))
    , gam(pf_in.gam)
    , mach(pf_in.mach)
    , reynolds(reynolds_in)
    , prandtl(prandtl_in)
    , heat_prefactor(1. / (reynolds * prandtl * (gam - 1) * mach * mach))
{
}

double DissipationTool::Txx(const CartesianGrid& grid, int ind) const
{
    return 1. / reynolds * (4 / 3. * der->DX(grid, alias::U, ind)
                               - 2 / 3. * der->DY(grid, alias::V, ind));
}

double DissipationTool::Txy(const CartesianGrid& grid, int ind) const
{
    return 1. / reynolds
        * (der->DX(grid, alias::V, ind) + der->DY(grid, alias::U, ind));
}

double DissipationTool::Tyy(const CartesianGrid& grid, int ind) const
{
    return 1. / reynolds * (4 / 3. * der->DY(grid, alias::V, ind)
                               - 2 / 3. * der->DX(grid, alias::U, ind));
}

double DissipationTool::dTxx_dx(const CartesianGrid& grid, int ind) const
{
    return 1. / reynolds * (4 / 3. * der->DXX(grid, alias::U, ind)
                               - 2 / 3. * der->DXY(grid, alias::V, ind));
}

double DissipationTool::dTxy_dx(const CartesianGrid& grid, int ind) const
{
    return 1. / reynolds
        * (der->DXY(grid, alias::U, ind) + der->DXX(grid, alias::V, ind));
}

double DissipationTool::dTxy_dy(const CartesianGrid& grid, int ind) const
{
    return 1. / reynolds
        * (der->DXY(grid, alias::V, ind) + der->DYY(grid, alias::U, ind));
}

double DissipationTool::dTyy_dy(const CartesianGrid& grid, int ind) const
{
    return 1. / reynolds * (4 / 3. * der->DYY(grid, alias::V, ind)
                               - 2 / 3. * der->DXY(grid, alias::U, ind));
}

double DissipationTool::dqx_dx(const CartesianGrid& grid, int ind) const
{
    return heat_prefactor * der->DXX(grid, alias::T, ind);
}

double DissipationTool::dqy_dy(const CartesianGrid& grid, int ind) const
{
    return heat_prefactor * der->DYY(grid, alias::T, ind);
}
