/**
 * \file irregular_derivatives.cpp
 * @brief Implementation of IrregularDerivatives class
 */
#include "irregular_derivatives.hpp"
#include "../grid/cartesian_grid.hpp"
#include "../grid/karagiozis_grid.hpp"
#include "../utils/flag_handler.hpp"
#include "../utils/flux_def.hpp"
#include "../utils/operators_overloads.hpp"

#include <iostream>
#include <utility>

template <typename Func>
auto IrregularDerivatives::zero_der() const
{
    return 0.0;
}

template <>
auto IrregularDerivatives::zero_der<PointToFlux>() const
{
    return Flux({0.0, 0.0, 0.0, 0.0});
}

template <>
auto IrregularDerivatives::zero_der<PositionToFlux>() const
{
    return Flux({0.0, 0.0, 0.0, 0.0});
}

template <typename Func>
auto IrregularDerivatives::first_der_two_discs(Func& func,
    GenericDiscontinuity* left, GenericDiscontinuity* right, double h) const
{
    const double eta = 1 - left->frac;
    const double eps = right->frac;
    if (eta + eps >= 1) { // Total space is larger than one grid space
        return (func(right->left()) - func(left->right())) / ((eps + eta) * h);
    }
    return zero_der<Func>();
}

template <typename Func>
auto IrregularDerivatives::first_der_one_disc(Func& func, const Point& q_disc,
    const Point& /*q1*/, const Point& q2, double eta, double h) const
{
    // Skewed-central derivative
    // Collapses to central when eta=1, and one-sided when eta=0
    return (func(q2) - func(q_disc)) / ((1. + eta) * h);
}

template <typename Func, GetDisc first_disc, GetDisc last_disc>
auto IrregularDerivatives::select_first_irregular_derivative(
    const KaragiozisGrid& grid, Func& func, int ind, int der_flag,
    HasDiscs& has_discs, int shift, double h) const
{
    if (has_discs.left and has_discs.right) { // Stuck between two discs
        auto disc_left = (grid.*last_disc)(ind - shift);
        auto disc_right = (grid.*first_disc)(ind);
        return first_der_two_discs<Func>(func, disc_left, disc_right, h);
    }
    if (has_discs.left
        and flag_functions::right(der_flag)
            != 0) { // Disc to the left but not a domain boundary to the right
        auto disc_left = (grid.*last_disc)(ind - shift);
        auto disc_point_left = disc_left->right();
        return first_der_one_disc<Func>(func, disc_point_left, grid.values(ind),
            grid.values(ind + shift), 1 - disc_left->frac, h);
    }
    if (has_discs.right
        and flag_functions::left(der_flag)
            != 0) { // Disc to the right but not a domain boundary to the left
        auto disc_right = (grid.*first_disc)(ind);
        auto disc_point_right = disc_right->left();
        return -first_der_one_disc<Func>(func, disc_point_right,
            grid.values(ind), grid.values(ind - shift), disc_right->frac, h);
    }
    return zero_der<Func>(); // Not enough space -> zero
}

IrregularDerivatives::IrregularDerivatives(PointFunctions& pf_in, int shiftX_in,
    int shiftY_in, double dx_in, double dy_in)
    : RegularDerivatives(pf_in, shiftX_in, shiftY_in, dx_in, dy_in)
{
}

auto IrregularDerivatives::disc_point_and_frac_left(
    GenericDiscontinuity* disc) const
{
    // Frac is always stored relative to grid point to the left
    // so the distance to the point at the right is 1-frac
    return std::make_tuple(disc->right(), 1 - disc->frac);
}

auto IrregularDerivatives::disc_point_and_frac_right(
    GenericDiscontinuity* disc) const
{
    return std::make_tuple(disc->left(), disc->frac);
}

template <typename Func>
auto IrregularDerivatives::DX_inner(
    const CartesianGrid& grid, Func func, int ind) const
{
    auto& in_grid = dynamic_cast<const KaragiozisGrid&>(grid);
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    bool has_left_disc = in_grid.has_discont_x(ind, -1);
    bool has_right_disc = in_grid.has_discont_x(ind);
    HasDiscs has_discs = {has_left_disc, has_right_disc};

    // Dispatch to irregular if have disc in imediate neighbors
    if (has_left_disc or has_right_disc) {
        return select_first_irregular_derivative<Func,
            &KaragiozisGrid::first_disc_x, &KaragiozisGrid::last_disc_x>(
            in_grid, func, ind, derx_flag, has_discs, shiftX, dx);
    }
    return RegularDerivatives::DX(grid, func, ind);
}

template <typename Func>
auto IrregularDerivatives::DY_inner(
    const CartesianGrid& grid, Func func, int ind) const
{
    auto& in_grid = dynamic_cast<const KaragiozisGrid&>(grid);
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    bool has_left_disc = in_grid.has_discont_y(ind, -1);
    bool has_right_disc = in_grid.has_discont_y(ind);
    HasDiscs has_discs = {has_left_disc, has_right_disc};

    // Dispatch to irregular if have disc in imediate neighbors
    if (has_left_disc or has_right_disc) {
        return select_first_irregular_derivative<Func,
            &KaragiozisGrid::first_disc_y, &KaragiozisGrid::last_disc_y>(
            in_grid, func, ind, dery_flag, has_discs, shiftY, dy);
    }
    return RegularDerivatives::DY(grid, func, ind);
}

double IrregularDerivatives::DX(
    const CartesianGrid& grid, alias::PointProperty func, int ind) const
{
    // Converts from PositionToDouble to PointToDouble
    auto f = [&](const Point& p) { return (pf.*func)(p); };
    return this->DX(grid, f, ind);
}

double IrregularDerivatives::DY(
    const CartesianGrid& grid, alias::PointProperty func, int ind) const
{
    // Converts from PositionToDouble to PointToDouble
    auto f = [&](const Point& p) { return (pf.*func)(p); };
    return this->DY(grid, f, ind);
}

double IrregularDerivatives::DX(
    const CartesianGrid& grid, PointToDouble& func, int ind) const
{
    return DX_inner<PointToDouble>(grid, func, ind);
}

double IrregularDerivatives::DY(
    const CartesianGrid& grid, PointToDouble& func, int ind) const
{
    return DY_inner<PointToDouble>(grid, func, ind);
}

Flux IrregularDerivatives::DX(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    return DX_inner<PointToFlux>(grid, func, ind);
}

Flux IrregularDerivatives::DXForward(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    auto& in_grid = dynamic_cast<const KaragiozisGrid&>(grid);
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    bool has_left_disc = in_grid.has_discont_x(ind, -1);
    bool has_right_disc = in_grid.has_discont_x(ind);
    // Because of numerical stability issues, if a discontinuity is found
    // in the neighborhood, cannot use split method (the only place where
    // Forward/Backward ders are used).
    // Either both fluxes are derived to one side, or none
    // So if possible will always dispatch to regular
    if (has_left_disc or has_right_disc
        or flag_functions::right(derx_flag) == 0) {
        return this->DX(grid, func, ind);
    }
    return RegularDerivatives::DXForward(grid, func, ind);
}

Flux IrregularDerivatives::DXBackward(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    auto& in_grid = dynamic_cast<const KaragiozisGrid&>(grid);
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    bool has_left_disc = in_grid.has_discont_x(ind, -1);
    bool has_right_disc = in_grid.has_discont_x(ind);
    // See DXForward
    if (has_left_disc or has_right_disc
        or flag_functions::left(derx_flag) == 0) {
        return this->DX(grid, func, ind);
    }
    return RegularDerivatives::DXBackward(grid, func, ind);
}

Flux IrregularDerivatives::DY(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    return DY_inner<PointToFlux>(grid, func, ind);
}

Flux IrregularDerivatives::DYForward(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    auto& in_grid = dynamic_cast<const KaragiozisGrid&>(grid);
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    bool has_left_disc = in_grid.has_discont_y(ind, -1);
    bool has_right_disc = in_grid.has_discont_y(ind);
    // See DXForward
    if (has_left_disc or has_right_disc
        or flag_functions::right(dery_flag) == 0) {
        return this->DY(grid, func, ind);
    }
    return RegularDerivatives::DYForward(grid, func, ind);
}

Flux IrregularDerivatives::DYBackward(
    const CartesianGrid& grid, PointToFlux& func, int ind) const
{
    auto& in_grid = dynamic_cast<const KaragiozisGrid&>(grid);
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    bool has_left_disc = in_grid.has_discont_y(ind, -1);
    bool has_right_disc = in_grid.has_discont_y(ind);
    // See DXForward
    if (has_left_disc or has_right_disc
        or flag_functions::left(dery_flag) == 0) {
        return this->DY(grid, func, ind);
    }
    return RegularDerivatives::DYBackward(grid, func, ind);
}

double IrregularDerivatives::DX(
    const CartesianGrid& grid, PositionToDouble& func, int ind) const
{
    // This is a function only used by dq_dx at Adiabatic Walls
    // It is NOT the correct way to compute the values
    return RegularDerivatives::DX(grid, func, ind);
}

double IrregularDerivatives::DY(
    const CartesianGrid& grid, PositionToDouble& func, int ind) const
{
    // This is a function only used by dq_dy at Adiabatic Walls
    // It is NOT the correct way to compute the values
    return RegularDerivatives::DY(grid, func, ind);
}

bool IrregularDerivatives::second_der_is_zero(
    HasDiscs& has_discs, int der_flag) const
{
    // Surrounded
    if (has_discs.left and has_discs.right) {
        return true;
    }
    // At left boundary with discont
    if (flag_functions::left(der_flag) == 0 and has_discs.right) {
        return true;
    }
    // At right boundary with discont
    if (flag_functions::right(der_flag) == 0 and has_discs.left) {
        return true;
    }
    // Width = 1 or 2
    if (flag_functions::left(der_flag) == 0
        and flag_functions::right(der_flag) < 2) {
        return true;
    }
    // Same as before, to be safe
    if (flag_functions::right(der_flag) == 0
        and flag_functions::left(der_flag) < 2) {
        return true;
    }
    return false;
}

int IrregularDerivatives::shift_pos_to_derive(int der_flag) const
{
    // Used in DXX and DYY
    // Second derivatives at the boundaries are the same as one point before. So
    // we check if we are at the boundary and shift to the next point to compute
    // the derivative there instead
    if (flag_functions::left(der_flag) == 0) { // Left boundary
        return +1;
    }
    if (flag_functions::right(der_flag) == 0) { // Right boundary
        return -1;
    }
    return 0; // Not a boundary
}

double IrregularDerivatives::DXX(
    const CartesianGrid& grid, alias::PointProperty func, int ind) const
{
    auto& in_grid = dynamic_cast<const KaragiozisGrid&>(grid);
    auto derx_flag = flag_functions::derx(grid.flag(ind));
    // Discs in immediate neighbourhood
    HasDiscs has_discs_close
        = {in_grid.has_discont_x(ind, -1), in_grid.has_discont_x(ind)};

    if (second_der_is_zero(has_discs_close, derx_flag)) {
        return 0.0;
    }
    int move_pos_to_derive = shiftX * shift_pos_to_derive(derx_flag);
    // Is on the border, so move one point to the side
    if (move_pos_to_derive != 0) {
        // Have to perform all checks again starting from new position
        return DXX(grid, func, ind + move_pos_to_derive);
    }
    // Is irregular
    if (has_discs_close.left or has_discs_close.right) {
        // Second derivative needs wider stencil, so check for discontinuities
        // two grid points away
        HasDiscs has_discs_far
            = {in_grid.has_discont_x(ind, -2), in_grid.has_discont_x(ind, 1)};
        // Convert to PointToDouble
        auto f = [&](const Point& p) { return (pf.*func)(p); };
        return select_second_irregular_derivative<&KaragiozisGrid::first_disc_x,
            &KaragiozisGrid::last_disc_x>(in_grid, f, ind, derx_flag,
            has_discs_close, has_discs_far, shiftX, dx);
    }
    return RegularDerivatives::DXX(grid, func, ind);
}

double IrregularDerivatives::DYY(
    const CartesianGrid& grid, alias::PointProperty func, int ind) const
{
    // Check DXX comments
    auto& in_grid = dynamic_cast<const KaragiozisGrid&>(grid);
    auto dery_flag = flag_functions::dery(grid.flag(ind));
    HasDiscs has_discs_close
        = {in_grid.has_discont_y(ind, -1), in_grid.has_discont_y(ind)};

    if (second_der_is_zero(has_discs_close, dery_flag)) {
        return 0.0;
    }
    int move_pos_to_derive = shiftY * shift_pos_to_derive(dery_flag);
    // Is on the border, so move one point to the side
    if (move_pos_to_derive != 0) {
        return DYY(grid, func, ind + move_pos_to_derive);
    }
    if (has_discs_close.left or has_discs_close.right) {
        HasDiscs has_discs_far
            = {in_grid.has_discont_y(ind, -2), in_grid.has_discont_y(ind, 1)};
        auto f = [&](const Point& p) { return (pf.*func)(p); };
        return select_second_irregular_derivative<&KaragiozisGrid::first_disc_y,
            &KaragiozisGrid::last_disc_y>(in_grid, f, ind, dery_flag,
            has_discs_close, has_discs_far, shiftY, dy);
    }
    return RegularDerivatives::DYY(grid, func, ind);
}

template <GetDisc first_disc, GetDisc last_disc>
double IrregularDerivatives::select_second_irregular_derivative(
    const KaragiozisGrid& grid, PointToDouble& func, int ind, int der_flag,
    HasDiscs& has_discs_close, HasDiscs& has_discs_far, int shift,
    double h) const
{
    // At this point we know for sure the point is irregular
    Point disc_point_left;
    Point disc_point_right;
    double eta; // Distance to disc on the left
    double eps; // Distance to disc on the right
    if (has_discs_close.left) {
        std::tie(disc_point_left, eta)
            = disc_point_and_frac_left((grid.*last_disc)(ind - shift));
        // Point before last at right boundary
        // Treat boundary point as if a discontinuity is on top of it
        if (flag_functions::right(der_flag) == 1) {
            return second_der_two_disc(func, disc_point_left, grid.values(ind),
                grid.values(ind + shift), grid.values(ind + shift), eta, 0.0,
                h);
        }
        if (has_discs_far.right) { // Actual two discs
            std::tie(disc_point_right, eps)
                = disc_point_and_frac_right((grid.*first_disc)(ind + shift));
            return second_der_two_disc(func, disc_point_left, grid.values(ind),
                grid.values(ind + shift), disc_point_right, eta, eps, h);
        }
        return second_der_one_disc(func, disc_point_left, grid.values(ind),
            grid.values(ind + shift), grid.values(ind + 2 * shift), eta, h);
    }
    if (has_discs_close.right) {
        std::tie(disc_point_right, eps)
            = disc_point_and_frac_right((grid.*first_disc)(ind));
        // Point before last at left boundary
        // Treat boundary point as if a discontinuity is on top of it
        if (flag_functions::left(der_flag) == 1) {
            return second_der_two_disc(func, grid.values(ind - shift),
                grid.values(ind - shift), grid.values(ind), disc_point_right,
                0.0, eps, h);
        }
        if (has_discs_far.left) { // Actual two discs
            std::tie(disc_point_left, eta)
                = disc_point_and_frac_left((grid.*last_disc)(ind - 2 * shift));
            return second_der_two_disc(func, disc_point_left,
                grid.values(ind - shift), grid.values(ind), disc_point_right,
                eta, eps, h);
        }
        return second_der_one_disc(func, disc_point_right, grid.values(ind),
            grid.values(ind - shift), grid.values(ind - 2 * shift), eps, h);
    }
    return 0.0;
}

double IrregularDerivatives::DXY(
    const CartesianGrid& grid, alias::PointProperty func, int ind) const
{
    auto& in_grid = dynamic_cast<const KaragiozisGrid&>(grid);
    HasDiscs has_discs_x
        = {in_grid.has_discont_x(ind, -1), in_grid.has_discont_x(ind)};
    HasDiscs has_discs_y
        = {in_grid.has_discont_y(ind, -1), in_grid.has_discont_y(ind)};

    int ind_left = grid.indJMinusOne(ind);
    HasDiscs left_point = {
        in_grid.has_discont_y(ind_left, -1), in_grid.has_discont_y(ind_left)};

    int ind_right = grid.indJPlusOne(ind);
    HasDiscs right_point = {
        in_grid.has_discont_y(ind_right, -1), in_grid.has_discont_y(ind_right)};

    int ind_bottom = grid.indIMinusOne(ind);
    HasDiscs bottom_point = {in_grid.has_discont_x(ind_bottom, -1),
        in_grid.has_discont_x(ind_bottom)};

    int ind_top = grid.indIPlusOne(ind);
    HasDiscs top_point
        = {in_grid.has_discont_x(ind_top, -1), in_grid.has_discont_x(ind_top)};

    /*
     * TL---TOP---TR
     * |     |     |
     * |     |     |
     * |     |     |
     * L----PNT----R
     * |     |     |
     * |     |     |
     * |     |     |
     * BL----B----BR
     */

    // Have to check if there is an open path from the point to the corner
    // If it has an open path, compute the estimate from cross_der_estimate()
    // Result is the mean of estimates divided by grid spacing

    int ind_corner;
    double der_estimate = 0;
    int counter = 0;
    // Bottom-left
    ind_corner = grid.IND(grid.indI(ind) - 1, grid.indJ(ind) - 1);
    if (ind_corner >= 0) {
        bool is_blocked = (has_discs_y.left or bottom_point.left)
            and (has_discs_x.left or left_point.left);
        if (not is_blocked) {
            der_estimate += cross_der_estimate(grid, func, ind, -1, -1);
            counter++;
        }
    }

    // Bottom-right
    ind_corner = grid.IND(grid.indI(ind) - 1, grid.indJ(ind) + 1);
    if (ind_corner >= 0) {
        bool is_blocked = (has_discs_y.left or bottom_point.right)
            and (has_discs_x.right or right_point.left);
        if (not is_blocked) {
            der_estimate += cross_der_estimate(grid, func, ind, -1, +1);
            counter++;
        }
    }

    // Top-left
    ind_corner = grid.IND(grid.indI(ind) + 1, grid.indJ(ind) - 1);
    if (ind_corner >= 0) {
        bool is_blocked = (has_discs_y.right or top_point.left)
            and (has_discs_x.left or left_point.right);
        if (not is_blocked) {
            der_estimate += cross_der_estimate(grid, func, ind, +1, -1);
            counter++;
        }
    }

    // Top-right
    ind_corner = grid.IND(grid.indI(ind) + 1, grid.indJ(ind) + 1);
    if (ind_corner >= 0) {
        bool is_blocked = (has_discs_y.right or top_point.right)
            and (has_discs_x.right or right_point.right);
        if (not is_blocked) {
            der_estimate += cross_der_estimate(grid, func, ind, +1, +1);
            counter++;
        }
    }
    if (counter == 0) {
        return 0.;
    }
    return der_estimate / (counter * dx * dy);
}

double IrregularDerivatives::second_der_one_disc(PointToDouble& func,
    const Point& q_disc, const Point& /*q1*/, const Point& q2, const Point& q3,
    double eta, double h) const
{
    // First order approximation from Taylor expansion. Always stable
    return 2 * (func(q_disc) - (eta + 2) * func(q2) + (eta + 1) * func(q3))
        / (h * h * (eta + 1) * (eta + 2));
}

double IrregularDerivatives::second_der_two_disc(PointToDouble& func,
    const Point& q_left, const Point& q1, const Point& q2, const Point& q_right,
    double eta, double eps, double h) const
{
    // Crafted second derivative approximation. Always underestimates the value
    // of the second derivative!!!!
    // O---x---O-------O-x-----O
    //            ^
    //         mid_val
    // mid_val is located at the mid-point of the discontinuities.
    // Its value is interpolated from the two grid points around it
    // Central second derivative is computed using the three values
    double mid_val = (eps - eta + 1) * (func(q2) - func(q1)) + 2 * func(q1);
    return 4 * (func(q_left) - mid_val + func(q_right))
        / (h * h * pow(eta + eps + 1, 2));
}

double IrregularDerivatives::cross_der_estimate(const CartesianGrid& grid,
    alias::PointProperty func, int ind, int i, int j) const
{
    // Estimate the cross derivative at the point ind from values at new point
    //(ind.I + i, ind.J + j) using Taylor expansion around the new point
    auto f = [&](const Point& p) { return (pf.*func)(p); };
    int corner = grid.IND(grid.indI(ind) + i, grid.indJ(ind) + j);
    double local_val = f(grid.values(ind));
    double corner_val = f(grid.values(corner));

    // This actualy returns dx*dy*cross_der
    return i * j * (corner_val - (local_val + j * dx * DX(grid, func, ind)
                                     + i * dy * DY(grid, func, ind)
                                     + dx * dx * DXX(grid, func, ind) / 2
                                     + dy * dy * DYY(grid, func, ind) / 2));
}
