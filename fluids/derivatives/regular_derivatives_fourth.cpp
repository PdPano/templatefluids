/**
 * \file regular_derivatives_fourth.cpp
 * @brief Implementation of RegularDerivativesFourth class
 */

#include "regular_derivatives_fourth.hpp"
#include "../grid/cartesian_grid.hpp"
#include "../utils/flag_handler.hpp"
#include "../utils/flux_def.hpp"
#include "../utils/operators_overloads.hpp"

#include <utility>

RegularDerivativesFourth::RegularDerivativesFourth(PointFunctions& pf_in,
    int shiftX_in, int shiftY_in, double dx_in, double dy_in)
    : Derivatives(pf_in, shiftX_in, shiftY_in, dx_in, dy_in)
{
}

template <typename Func>
auto RegularDerivativesFourth::first_derivative_central(
    Func& func, int ind, int shift, double h) const
{
    return (func(ind - 2 * shift) - 8 * func(ind - shift)
               + 8 * func(ind + shift) - func(ind + 2 * shift))
        / (12 * h);
}

template <typename Func>
auto RegularDerivativesFourth::first_derivative_left(
    Func& func, int ind, int shift, double h) const
{
    return first_derivative_right<Func>(func, ind, -shift, -h);
}

template <typename Func>
auto RegularDerivativesFourth::first_derivative_right(
    Func& func, int ind, int shift, double h) const
{
    return (-11 * func(ind) + 18 * func(ind + shift) - 9 * func(ind + 2 * shift)
               + 2 * func(ind + 3 * shift))
        / (6 * h);
}

template <typename Func>
auto RegularDerivativesFourth::first_derivative_almost_left(
    Func& func, int ind, int shift, double h) const
{
    return first_derivative_almost_right<Func>(func, ind, -shift, -h);
}

template <typename Func>
auto RegularDerivativesFourth::first_derivative_almost_right(
    Func& func, int ind, int shift, double h) const
{
    return (-3 * func(ind - shift) - 10 * func(ind) + 18 * func(ind + shift)
               - 6 * func(ind + 2 * shift) + func(ind + 3 * shift))
        / (12 * h);
}

template <typename Func>
auto RegularDerivativesFourth::zero_der() const
{
    return 0.0;
}

template <>
auto RegularDerivativesFourth::zero_der<PointToFlux>() const
{
    return Flux({0.0, 0.0, 0.0, 0.0});
}

template <>
auto RegularDerivativesFourth::zero_der<PositionToFlux>() const
{
    return Flux({0.0, 0.0, 0.0, 0.0});
}
template <typename Func>
auto RegularDerivativesFourth::first_der_selector(
    Func& func, int ind, int shift, double h, uint16_t der_flag) const
{
    // Checks if the point has enough space on either side to calculate the
    // derivative
    if (std::min(
            flag_functions::left(der_flag), flag_functions::right(der_flag))
        > 1) {
        return first_derivative_central<Func>(func, ind, shift, h);
    }
    if (flag_functions::left(der_flag) == 0
        and flag_functions::right(der_flag) > 3) { // Left domain boundary
        return first_derivative_right<Func>(func, ind, shift, h);
    }
    if (flag_functions::left(der_flag) > 3
        and flag_functions::right(der_flag) == 0) { // Right domain boundary
        return first_derivative_left<Func>(func, ind, shift, h);
    }
    if (flag_functions::left(der_flag) == 1
        and flag_functions::right(der_flag) > 2) { // Left domain boundary
        return first_derivative_almost_right<Func>(func, ind, shift, h);
    }
    if (flag_functions::left(der_flag) > 2
        and flag_functions::right(der_flag) == 1) { // Right domain boundary
        return first_derivative_almost_left<Func>(func, ind, shift, h);
    }
    return zero_der<Func>();
}

double RegularDerivativesFourth::DX(
    const CartesianGrid& grid, alias::PointProperty func, int ind) const
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    // Converts from PointToDouble to PositionToDouble
    auto f = [&](int ind) { return (pf.*func)(grid.values(ind)); };
    return first_der_selector<PositionToDouble>(f, ind, shiftX, dx, derx_flag);
}

double RegularDerivativesFourth::DY(
    const CartesianGrid& grid, alias::PointProperty func, int ind) const
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    // Converts from PointToDouble to PositionToDouble
    auto f = [&](int ind) { return (pf.*func)(grid.values(ind)); };
    return first_der_selector<PositionToDouble>(f, ind, shiftY, dy, dery_flag);
}

double RegularDerivativesFourth::DX(
    const CartesianGrid& grid, PointToDouble& func, int ind) const
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    // Converts from PointToDouble to PositionToDouble
    auto f = [&](int ind) { return func(grid.values(ind)); };
    return first_der_selector<PositionToDouble>(f, ind, shiftX, dx, derx_flag);
}

double RegularDerivativesFourth::DY(
    const CartesianGrid& grid, PointToDouble& func, int ind) const
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    // Converts from PointToDouble to PositionToDouble
    auto f = [&](int ind) { return func(grid.values(ind)); };
    return first_der_selector<PositionToDouble>(f, ind, shiftY, dy, dery_flag);
}

double RegularDerivativesFourth::DX(
    const CartesianGrid& grid, PositionToDouble& func, int ind) const
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    return first_der_selector<PositionToDouble>(
        func, ind, shiftX, dx, derx_flag);
}
double RegularDerivativesFourth::DY(
    const CartesianGrid& grid, PositionToDouble& func, int ind) const
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    return first_der_selector<PositionToDouble>(
        func, ind, shiftY, dy, dery_flag);
}

Flux RegularDerivativesFourth::DX(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    // Converts from PointToFlux to PositionToFlux
    auto f = [&](int i) { return func(grid.values(i)); };
    return first_der_selector<PositionToFlux>(f, ind, shiftX, dx, derx_flag);
}

Flux RegularDerivativesFourth::DXForward(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    if (std::min(
            flag_functions::left(derx_flag), flag_functions::right(derx_flag))
        > 3) { // Is possible to compute
        // Converts from PointToFlux to PositionToFlux
        auto f = [&](int i) { return func(grid.values(i)); };
        return first_derivative_almost_right<PositionToFlux>(
            f, ind, shiftX, dx);
    }
    return this->DX(grid, func, ind);
}

Flux RegularDerivativesFourth::DXBackward(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    if (std::min(
            flag_functions::left(derx_flag), flag_functions::right(derx_flag))
        > 3) { // Is possible to compute
        // Converts from PointToFlux to PositionToFlux
        auto f = [&](int i) { return func(grid.values(i)); };
        return first_derivative_almost_left<PositionToFlux>(f, ind, shiftX, dx);
    }
    return this->DX(grid, func, ind);
}

Flux RegularDerivativesFourth::DY(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    // Converts from PointToFlux to PositionToFlux
    auto f = [&](int i) { return func(grid.values(i)); };
    return first_der_selector<PositionToFlux>(f, ind, shiftY, dy, dery_flag);
}

Flux RegularDerivativesFourth::DYForward(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    if (std::min(
            flag_functions::left(dery_flag), flag_functions::right(dery_flag))
        > 3) { // Is possible to compute
        // Converts from PointToFlux to PositionToFlux
        auto f = [&](int i) { return func(grid.values(i)); };
        return first_derivative_almost_right<PositionToFlux>(
            f, ind, shiftY, dy);
    }
    return this->DY(grid, func, ind);
}

Flux RegularDerivativesFourth::DYBackward(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    if (std::min(
            flag_functions::left(dery_flag), flag_functions::right(dery_flag))
        > 3) { // Is possible to compute
        // Converts from PointToFlux to PositionToFlux
        auto f = [&](int i) { return func(grid.values(i)); };
        return first_derivative_almost_left<PositionToFlux>(f, ind, shiftY, dy);
    }
    return this->DY(grid, func, ind);
}

Flux RegularDerivativesFourth::DX(
    const CartesianGrid& grid, PositionToFlux& func, int ind) const
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    return first_der_selector<PositionToFlux>(func, ind, shiftX, dx, derx_flag);
}

Flux RegularDerivativesFourth::DXForward(
    const CartesianGrid& grid, PositionToFlux& func, int ind) const
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    if (std::min(
            flag_functions::left(derx_flag), flag_functions::right(derx_flag))
        > 3) { // Is possible to compute
        return first_derivative_almost_right<PositionToFlux>(
            func, ind, shiftX, dx);
    }
    return this->DX(grid, func, ind);
}

Flux RegularDerivativesFourth::DXBackward(
    const CartesianGrid& grid, PositionToFlux& func, int ind) const
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    if (std::min(
            flag_functions::left(derx_flag), flag_functions::right(derx_flag))
        > 3) { // Is possible to compute
        return first_derivative_almost_left<PositionToFlux>(
            func, ind, shiftX, dx);
    }
    return this->DX(grid, func, ind);
}

Flux RegularDerivativesFourth::DY(
    const CartesianGrid& grid, PositionToFlux& func, int ind) const
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    return first_der_selector<PositionToFlux>(func, ind, shiftY, dy, dery_flag);
}

Flux RegularDerivativesFourth::DYForward(
    const CartesianGrid& grid, PositionToFlux& func, int ind) const
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    if (std::min(
            flag_functions::left(dery_flag), flag_functions::right(dery_flag))
        > 3) { // Is possible to compute
        return first_derivative_almost_right<PositionToFlux>(
            func, ind, shiftY, dy);
    }
    return this->DY(grid, func, ind);
}

Flux RegularDerivativesFourth::DYBackward(
    const CartesianGrid& grid, PositionToFlux& func, int ind) const
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    if (std::min(
            flag_functions::left(dery_flag), flag_functions::right(dery_flag))
        > 3) { // Is possible to compute
        return first_derivative_almost_left<PositionToFlux>(
            func, ind, shiftY, dy);
    }
    return this->DY(grid, func, ind);
}

double RegularDerivativesFourth::DXX(
    const CartesianGrid& grid, alias::PointProperty func, int ind) const
{
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    // Converts from PointToDouble to PositionToDouble
    auto f = [&](int ind) { return (pf.*func)(grid.values(ind)); };
    return second_der_selector(f, ind, shiftX, dx, derx_flag);
}

double RegularDerivativesFourth::DYY(
    const CartesianGrid& grid, alias::PointProperty func, int ind) const
{
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    // Converts from PointToDouble to PositionToDouble
    auto f = [&](int ind) { return (pf.*func)(grid.values(ind)); };
    return second_der_selector(f, ind, shiftY, dy, dery_flag);
}

double RegularDerivativesFourth::DXY(
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

double RegularDerivativesFourth::second_der_selector(PositionToDouble& func,
    int ind, int shift, double h, uint16_t der_flag) const
{
    // Checks if the point has enough space on either side to calculate the
    // derivative
    if (std::min(
            flag_functions::left(der_flag), flag_functions::right(der_flag))
        > 1) {
        return second_derivative_central(func, ind, shift, h);
    }
    if (flag_functions::left(der_flag) == 0
        and flag_functions::right(der_flag) >= 4) {
        return second_derivative_right(func, ind, shift, h);
    }
    if (flag_functions::left(der_flag) >= 4
        and flag_functions::right(der_flag) == 0) {
        return second_derivative_left(func, ind, shift, h);
    }
    if (flag_functions::left(der_flag) == 1
        and flag_functions::right(der_flag) >= 3) {
        return second_derivative_almost_right(func, ind, shift, h);
    }
    if (flag_functions::left(der_flag) >= 3
        and flag_functions::right(der_flag) == 1) {
        return second_derivative_almost_left(func, ind, shift, h);
    }
    return 0.0;
}

double RegularDerivativesFourth::second_derivative_central(
    PositionToDouble& func, int ind, int shift, double h) const
{
    return (-func(ind - 2 * shift) + 16 * func(ind - shift) - 30 * func(ind)
               + 16 * func(ind + shift) - func(ind + 2 * shift))
        / (12 * h * h);
}

double RegularDerivativesFourth::second_derivative_right(
    PositionToDouble& func, int ind, int shift, double h) const
{
    return (35 / 12. * func(ind) - 26 / 3. * func(ind + shift)
               + 19 / 2. * func(ind + 2 * shift)
               + -14 / 3. * func(ind + 3 * shift)
               + 11 / 12. * func(ind + 4 * shift))
        / (h * h);
}

double RegularDerivativesFourth::second_derivative_left(
    PositionToDouble& func, int ind, int shift, double h) const
{
    return second_derivative_right(func, ind, -shift, h);
}

double RegularDerivativesFourth::second_derivative_almost_left(
    PositionToDouble& func, int ind, int shift, double h) const
{
    return second_derivative_almost_right(func, ind, -shift, h);
}

double RegularDerivativesFourth::second_derivative_almost_right(
    PositionToDouble& func, int ind, int shift, double h) const
{
    return (10 * func(ind - shift) - 15 * func(ind) - 4 * func(ind + shift)
               + 14 * func(ind + 2 * shift) - 6 * func(ind + 3 * shift)
               + func(ind + 4 * shift))
        / (12 * h * h);
}
