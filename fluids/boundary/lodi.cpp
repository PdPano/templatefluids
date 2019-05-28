/**
 * \file lodi.cpp
 * @brief Implementation of Lodi class
 */
#include "lodi.hpp"
#include "../grid/cartesian_grid.hpp"
#include "../utils/point_def.hpp"
#include "../utils/point_functions.hpp"

#include <utility>
Lodi::Lodi(PointFunctions& pf_in, std::shared_ptr<Derivatives> der_in)
    : pf(pf_in)
    , der(std::move(der_in))
{
}

LodiArray Lodi::values(const CartesianGrid& grid, DerFunction df,
    alias::PointProperty vel_1, alias::PointProperty vel_2, int ind) const
{
    LodiArray lodi{}; // std::array<double, 4>

    auto point = grid.values(ind);
    double u = (pf.*vel_1)(point);
    double c = pf.sound_speed(point);
    double rho = pf.rho(point);

    auto Der = [&](alias::PointProperty func) {
        return ((der.get())->*df)(grid, func, ind);
    };

    lodi[0] = (u - c) * (Der(alias::P) - rho * c * Der(vel_1));
    lodi[1] = (u) * (c * c * Der(alias::RHO) - Der(alias::P));
    lodi[2] = (u) * (Der(vel_2));
    lodi[3] = (u + c) * (Der(alias::P) + rho * c * Der(vel_1));

    return lodi;
}
