/**
 * \file regular_derivatives.cpp
 * @brief Implementation of RegularDerivatives class
 */

#include "regular_derivatives.hpp"
#include "../grid/cartesian_grid.hpp"
#include "../utils/flag_handler.hpp"
#include "../utils/flux_def.hpp"
#include "../utils/operators_overloads.hpp"

#include <utility>

RegularDerivatives::RegularDerivatives(PointFunctions& pf_in, int shiftX_in,
    int shiftY_in, double dx_in, double dy_in)
    : Derivatives(pf_in, shiftX_in, shiftY_in, dx_in, dy_in)
{
}

template <typename Func>
auto RegularDerivatives::first_derivative_central(
    Func& func, int ind, int shift, double h) const
{
    return (func(ind + shift) - func(ind - shift)) / (2 * h);
}

template <typename Func>
auto RegularDerivatives::first_derivative_right(
    Func& func, int ind, int shift, double h) const
{
    return (func(ind + shift) - func(ind)) / h;
}

template <typename Func>
auto RegularDerivatives::first_derivative_left(
    Func& func, int ind, int shift, double h) const
{
    // Equivalent to (func(ind) - func(ind-shift))/h
    return first_derivative_right<Func>(func, ind, -shift, -h);
}

template <typename Func>
auto RegularDerivatives::zero_der() const
{
    return 0.0;
}

template <>
auto RegularDerivatives::zero_der<PointToFlux>() const
{
    return Flux({0.0, 0.0, 0.0, 0.0});
}

template <>
auto RegularDerivatives::zero_der<PositionToFlux>() const
{
    return Flux({0.0, 0.0, 0.0, 0.0});
}

template <typename Func>
auto RegularDerivatives::first_der_selector(
    Func& func, int ind, int shift, double h, int der_flag) const
{
    // Checks if the point has enough space on either side to calculate the
    // derivative
    if (std::min(
            flag_functions::left(der_flag), flag_functions::right(der_flag))
        != 0) { // Both are >0
        return first_derivative_central<Func>(func, ind, shift, h);
    }
    if (flag_functions::left(der_flag) == 0
        and flag_functions::right(der_flag) != 0) { // Left domain boundary
        return first_derivative_right<Func>(func, ind, shift, h);
    }
    if (flag_functions::left(der_flag) != 0
        and flag_functions::right(der_flag) == 0) { // Right domain boundary
        return first_derivative_left<Func>(func, ind, shift, h);
    }
    return zero_der<Func>();
}

double RegularDerivatives::DX(
    const CartesianGrid& grid, alias::PointProperty func, int ind) const
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    // Converts from PointToDouble to PositionToDouble
    auto f = [&](int ind) { return (pf.*func)(grid.values(ind)); };
    return first_der_selector<PositionToDouble>(f, ind, shiftX, dx, derx_flag);
}

double RegularDerivatives::DY(
    const CartesianGrid& grid, alias::PointProperty func, int ind) const
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    // Converts from PointToDouble to PositionToDouble
    auto f = [&](int ind) { return (pf.*func)(grid.values(ind)); };
    return first_der_selector<PositionToDouble>(f, ind, shiftY, dy, dery_flag);
}

double RegularDerivatives::DX(
    const CartesianGrid& grid, PointToDouble& func, int ind) const
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    // Converts from PointToDouble to PositionToDouble
    auto f = [&](int ind) { return func(grid.values(ind)); };
    return first_der_selector<PositionToDouble>(f, ind, shiftX, dx, derx_flag);
}

double RegularDerivatives::DY(
    const CartesianGrid& grid, PointToDouble& func, int ind) const
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    // Converts from PointToDouble to PositionToDouble
    auto f = [&](int ind) { return func(grid.values(ind)); };
    return first_der_selector<PositionToDouble>(f, ind, shiftY, dy, dery_flag);
}

double RegularDerivatives::DX(
    const CartesianGrid& grid, PositionToDouble& func, int ind) const
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    return first_der_selector<PositionToDouble>(
        func, ind, shiftX, dx, derx_flag);
}
double RegularDerivatives::DY(
    const CartesianGrid& grid, PositionToDouble& func, int ind) const
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    return first_der_selector<PositionToDouble>(
        func, ind, shiftY, dy, dery_flag);
}

Flux RegularDerivatives::DX(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    // Converts from PointToFlux to PositionToFlux
    auto f = [&](int i) { return func(grid.values(i)); };
    return first_der_selector<PositionToFlux>(f, ind, shiftX, dx, derx_flag);
}

Flux RegularDerivatives::DXForward(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    if (flag_functions::right(derx_flag) != 0) { // Is possible to compute
        // Converts from PointToFlux to PositionToFlux
        auto f = [&](int i) { return func(grid.values(i)); };
        return first_derivative_right<PositionToFlux>(f, ind, shiftX, dx);
    }
    return this->DX(grid, func, ind);
}

Flux RegularDerivatives::DXBackward(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    if (flag_functions::left(derx_flag) != 0) { // Is possible to compute
        // Converts from PointToFlux to PositionToFlux
        auto f = [&](int i) { return func(grid.values(i)); };
        return first_derivative_left<PositionToFlux>(f, ind, shiftX, dx);
    }
    return this->DX(grid, func, ind);
}

Flux RegularDerivatives::DY(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    // Converts from PointToFlux to PositionToFlux
    auto f = [&](int i) { return func(grid.values(i)); };
    return first_der_selector<PositionToFlux>(f, ind, shiftY, dy, dery_flag);
}

Flux RegularDerivatives::DYForward(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    if (flag_functions::right(dery_flag) != 0) { // Is possible to compute
        // Converts from PointToFlux to PositionToFlux
        auto f = [&](int i) { return func(grid.values(i)); };
        return first_derivative_right<PositionToFlux>(f, ind, shiftY, dy);
    }
    return this->DY(grid, func, ind);
}

Flux RegularDerivatives::DYBackward(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    if (flag_functions::left(dery_flag) != 0) { // Is possible to compute
        // Converts from PointToFlux to PositionToFlux
        auto f = [&](int i) { return func(grid.values(i)); };
        return first_derivative_left<PositionToFlux>(f, ind, shiftY, dy);
    }
    return this->DY(grid, func, ind);
}

Flux RegularDerivatives::DX(const CartesianGrid& grid, PositionToFlux& func,
    int ind) const //!<First derivative in the X direction
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    return first_der_selector<PositionToFlux>(func, ind, shiftX, dx, derx_flag);
}

Flux RegularDerivatives::DXForward(const CartesianGrid& grid,
    PositionToFlux& func,
    int ind) const //!< Right-sided first derivative in the X direction
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    if (flag_functions::right(derx_flag) != 0) { // Is possible to compute
        return first_derivative_right<PositionToFlux>(func, ind, shiftX, dx);
    }
    return this->DX(grid, func, ind);
}

Flux RegularDerivatives::DXBackward(const CartesianGrid& grid,
    PositionToFlux& func,
    int ind) const //!< Left-sided first derivative in the X direction
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    if (flag_functions::left(derx_flag) != 0) { // Is possible to compute
        return first_derivative_left<PositionToFlux>(func, ind, shiftX, dx);
    }
    return this->DX(grid, func, ind);
}

Flux RegularDerivatives::DY(const CartesianGrid& grid, PositionToFlux& func,
    int ind) const //!<First derivative in the Y direction
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    return first_der_selector<PositionToFlux>(func, ind, shiftY, dy, dery_flag);
}

Flux RegularDerivatives::DYForward(const CartesianGrid& grid,
    PositionToFlux& func,
    int ind) const //!<Right-sided first derivative in the Y direction
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    if (flag_functions::right(dery_flag) != 0) { // Is possible to compute
        return first_derivative_right<PositionToFlux>(func, ind, shiftY, dy);
    }
    return this->DY(grid, func, ind);
}

Flux RegularDerivatives::DYBackward(const CartesianGrid& grid,
    PositionToFlux& func,
    int ind) const //!<Left-sided first derivative in the Y direction
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    if (flag_functions::left(dery_flag) != 0) { // Is possible to compute
        return first_derivative_left<PositionToFlux>(func, ind, shiftY, dy);
    }
    return this->DY(grid, func, ind);
}

double RegularDerivatives::DXX(
    const CartesianGrid& grid, alias::PointProperty func, int ind) const
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    // Converts from PointToDouble to PositionToDouble
    auto f = [&](int ind) { return (pf.*func)(grid.values(ind)); };
    return second_der_selector(f, ind, shiftX, dx, derx_flag);
}

double RegularDerivatives::DYY(
    const CartesianGrid& grid, alias::PointProperty func, int ind) const
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    // Converts from PointToDouble to PositionToDouble
    auto f = [&](int ind) { return (pf.*func)(grid.values(ind)); };
    return second_der_selector(f, ind, shiftY, dy, dery_flag);
}

double RegularDerivatives::DXY(
    const CartesianGrid& grid, alias::PointProperty func, int ind) const
{
    auto flag = grid.flag(ind);
    auto derx_flag = flag_functions::derx(flag);
    auto dery_flag = flag_functions::dery(flag);
    // Converts from PointToDouble to PositionToDouble
    auto f = [&](int ind) { return (pf.*func)(grid.values(ind)); };
    // Creates cross derivative from two first order derivatives
    auto f_x = [&](
        int ind) { return first_der_selector(f, ind, shiftX, dx, derx_flag); };
    auto f_y = [&](
        int ind) { return first_der_selector(f, ind, shiftY, dy, dery_flag); };

    // computes 1/2*(d^2 f / dx dy + d^2 f /dy dx)
    // This is necessary to avoid problems near inner corners as those in a step
    double res = (first_der_selector(f_x, ind, shiftY, dy, dery_flag)
                     + first_der_selector(f_y, ind, shiftX, dx, derx_flag))
        / 2;
    return res;
}

double RegularDerivatives::second_der_selector(
    PositionToDouble& func, int ind, int shift, double h, int der_flag) const
{
    // Checks if the point has enough space on either side to calculate the
    // derivative
    if (std::min(
            flag_functions::left(der_flag), flag_functions::right(der_flag))
        != 0) {
        return second_derivative_central(func, ind, shift, h);
    }
    if (flag_functions::left(der_flag) == 0
        and flag_functions::right(der_flag) >= 2) {
        return second_derivative_right(func, ind, shift, h);
    }
    if (flag_functions::left(der_flag) >= 2
        and flag_functions::right(der_flag) == 0) {
        return second_derivative_left(func, ind, shift, h);
    }
    return 0.0;
}

double RegularDerivatives::second_derivative_central(
    PositionToDouble& func, int ind, int shift, double h) const
{
    return (func(ind - shift) - 2 * func(ind) + func(ind + shift)) / (h * h);
}

double RegularDerivatives::second_derivative_right(
    PositionToDouble& func, int ind, int shift, double h) const
{
    return second_derivative_central(func, ind + shift, shift, h);
}

double RegularDerivatives::second_derivative_left(
    PositionToDouble& func, int ind, int shift, double h) const
{
    return second_derivative_central(func, ind - shift, shift, h);
}
