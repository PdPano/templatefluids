/**
 * \file boundary.cpp
 * @brief Implementation of Boundary class
 */
#include "boundary.hpp"
#include "../derivatives/derivatives.hpp"
#include "../grid/cartesian_grid.hpp"
#include "../reconstructions/dissipation_tool.hpp"
#include "../utils/boundary_point_def.hpp"
#include "../utils/global_vars.hpp"
#include "../utils/point_functions.hpp"

Boundary::Boundary(PointFunctions& pf_in, std::shared_ptr<Derivatives> der_in,
    double reynolds_in, double prandtl_in)
    : pf(pf_in)
    , der(der_in)
    , gam(pf_in.gam)
    , lodi(Lodi(pf_in, der_in))
    , dissipation_tool(std::make_shared<DissipationTool>(
          pf_in, der_in, reynolds_in, prandtl_in))
    , heat_prefactor(1. / (reynolds_in * prandtl_in * pf_in.mach * pf_in.mach
                              * (pf_in.gam - 1)))
{
}

Boundary::BoundaryConfiguration Boundary::configure(
    char direction, const BoundaryPoint& bp, const double& t) const
{
    BoundaryConfiguration bconf(bp);

    if (direction == 'X') {
        bconf.is_left_or_bottom = bp.is_left();
        bconf.v1 = alias::U;
        bconf.v2 = alias::V;
        bconf.df = &Derivatives::DX;
        bconf.t = t;
    }
    else if (direction == 'Y') {
        bconf.is_left_or_bottom = bp.is_bottom();
        bconf.v1 = alias::V;
        bconf.v2 = alias::U;
        bconf.df = &Derivatives::DY;
        bconf.t = t;
    }
    return bconf;
}

void Boundary::fix_boundary(
    CartesianGrid* grid, const BoundaryPoint& bp, const double& t) const
{
    Point p_values = grid->values(bp.ind); // Temporary value
    double c = pf.sound_speed(p_values);
    double u = fabs(pf.u(p_values));
    double v = fabs(pf.v(p_values));

    if (bp.is_subsonic_inlet()) {
        // If it is supersonic in the direction perpendicular to the boundary
        // it must also set the density
        if ((bp.is_x() and u > c) or (bp.is_y() and v > c)) {
            p_values.set_rho(bp.rho);
        }
        p_values.set_ru(p_values.rho() * bp.u * bp.time_function(t));
        p_values.set_rv(p_values.rho() * bp.v * bp.time_function(t));
        p_values.set_e(pf.e_from_T(p_values, bp.T));
    }

    if (bp.is_supersonic_inlet()) {
        // All values are fixed to the ones prescribed by the boundary
        // configuration
        p_values.set_rho(bp.rho);
        p_values.set_ru(p_values.rho() * bp.u * bp.time_function(t));
        p_values.set_rv(p_values.rho() * bp.v * bp.time_function(t));
        p_values.set_e(pf.e_from_T(p_values, bp.T));
    }

    // Correction comes from lodi relations only. It's a soft B.C.
    if (bp.is_subsonic_outlet()) {
    }

    // No corrections needed. Everything comes from within the domain
    if (bp.is_supersonic_outlet()) {
    }

    // Walls are assumed static
    if (bp.is_adiabatic_no_slip_wall()) {
        p_values.set_ru(0.);
        p_values.set_rv(0.);
    }

    if (bp.is_isotermal_no_slip_wall()) {
        p_values.set_ru(0.);
        p_values.set_rv(0.);
        p_values.set_e(pf.e_from_T(p_values, bp.T));
    }

    grid->set_values(p_values, bp.ind); // Update grid values
}

Flux Boundary::convection_x(
    const CartesianGrid& grid, const BoundaryPoint& bp, const double& t) const
{
    // Obtain the lodi values. Cannot use BoundaryConfiguration because
    // Lodi is out of scope.
    auto lodi_a
        = lodi.values(grid, &Derivatives::DX, alias::U, alias::V, bp.ind);
    auto bconf = configure('X', bp, t);

    // Select correction applied to lodi
    if (bp.is_subsonic_inlet_x()) {
        subsonic_inlet(&lodi_a, grid, bconf);
    }
    else if (bp.is_supersonic_inlet_x()) {
        supersonic_inlet(&lodi_a, grid, bconf);
    }
    else if (bp.is_subsonic_outlet_x()) {
        subsonic_outlet(&lodi_a, grid, bconf);
    }
    else if (bp.is_supersonic_outlet_x()) {
        supersonic_outlet(&lodi_a, grid, bconf);
    }
    else if (bp.is_adiabatic_no_slip_wall_x()) {
        adiabatic_no_slip_wall(&lodi_a, grid, bconf);
    }
    else if (bp.is_isotermal_no_slip_wall_x()) {
        isotermal_no_slip_wall(&lodi_a, grid, bconf);
    }

    // Compute flux from corrected lodi values
    return flux_from_lodi(grid, lodi_a, bconf);
}

Flux Boundary::convection_y(
    const CartesianGrid& grid, const BoundaryPoint& bp, const double& t) const
{
    // Same as convection_x
    auto lodi_a
        = lodi.values(grid, &Derivatives::DY, alias::V, alias::U, bp.ind);
    auto bconf = configure('Y', bp, t);

    if (bp.is_subsonic_inlet_y()) {
        subsonic_inlet(&lodi_a, grid, bconf);
    }
    if (bp.is_supersonic_inlet_y()) {
        supersonic_inlet(&lodi_a, grid, bconf);
    }
    if (bp.is_subsonic_outlet_y()) {
        subsonic_outlet(&lodi_a, grid, bconf);
    }
    if (bp.is_supersonic_outlet_y()) {
        supersonic_outlet(&lodi_a, grid, bconf);
    }
    if (bp.is_adiabatic_no_slip_wall_y()) {
        adiabatic_no_slip_wall(&lodi_a, grid, bconf);
    }
    if (bp.is_isotermal_no_slip_wall_y()) {
        isotermal_no_slip_wall(&lodi_a, grid, bconf);
    }

    return rotate_to_y(flux_from_lodi(grid, lodi_a, bconf));
}

void Boundary::subsonic_inlet(LodiArray* lodi_a, const CartesianGrid& grid,
    const BoundaryConfiguration& bconf) const
{
    auto point = grid.values(bconf.bp.ind);
    double c = pf.sound_speed(point);
    double rho = pf.rho(point);
    double time_der = bconf.bp.der_time_function(bconf.t);
    (*lodi_a)[2] = -((pf.*(bconf.v2))(point)) * time_der;
    if (bconf.is_left_or_bottom) {
        (*lodi_a)[3]
            = (*lodi_a)[0] - 2 * rho * c * ((pf.*(bconf.v1))(point)) * time_der;
    }
    else {
        (*lodi_a)[0]
            = (*lodi_a)[3] + 2 * rho * c * ((pf.*(bconf.v1))(point)) * time_der;
    }
}

void Boundary::supersonic_inlet(LodiArray* lodi_a,
    const CartesianGrid& /*grid*/, const BoundaryConfiguration& /*bconf*/) const
{
    // No flux -> all values fixed
    (*lodi_a) = {{0., 0., 0., 0.}};
}

void Boundary::subsonic_outlet(LodiArray* lodi_a, const CartesianGrid& grid,
    const BoundaryConfiguration& bconf) const
{
    auto point = grid.values(bconf.bp.ind);
    double vel_1 = (pf.*(bconf.v1))(point);
    double c = pf.sound_speed(point);

    if (fabs(vel_1) < c) { // Is in fact subsonic
        double pressure_correction = c * 0.9
            * (1 - max_mach_number * max_mach_number)
            * (pf.pressure(point) - bconf.bp.p);
        if (bconf.is_left_or_bottom) {
            (*lodi_a)[3] = pressure_correction;
        }
        else {
            (*lodi_a)[0] = pressure_correction;
        }
    }
}
void Boundary::supersonic_outlet(LodiArray* /*lodi_a*/,
    const CartesianGrid& /*grid*/, const BoundaryConfiguration& /*bconf*/) const
{
    // No changes -> everything computed from inside
}

void Boundary::adiabatic_no_slip_wall(LodiArray* lodi_a,
    const CartesianGrid& /*grid*/, const BoundaryConfiguration& bconf) const
{
    (*lodi_a)[1] = 0.0;
    (*lodi_a)[2] = 0.0;

    if (bconf.is_left_or_bottom) {
        (*lodi_a)[3] = (*lodi_a)[0];
    }
    else {
        (*lodi_a)[0] = (*lodi_a)[3];
    }
}

void Boundary::isotermal_no_slip_wall(LodiArray* lodi_a,
    const CartesianGrid& grid, const BoundaryConfiguration& bconf) const
{
    auto point = grid.values(bconf.bp.ind);
    double time_der = bconf.bp.der_time_function(bconf.t);
    (*lodi_a)[1] = 0.0;
    (*lodi_a)[2] = -((pf.*(bconf.v2))(point)) * time_der;

    if (bconf.is_left_or_bottom) {
        (*lodi_a)[3] = (*lodi_a)[0];
    }
    else {
        (*lodi_a)[0] = (*lodi_a)[3];
    }
}

Flux Boundary::dissipation_x(
    const CartesianGrid& grid, const BoundaryPoint& bp) const
{
    int ind = bp.ind;
    double Txx = dissipation_tool->Txx(grid, ind);
    double Txy = dissipation_tool->Txy(grid, ind);
    double dTxx_dx = dissipation_tool->dTxx_dx(grid, ind);
    double dTxy_dx = dissipation_tool->dTxy_dx(grid, ind);
    double dqx_dx = dissipation_tool->dqx_dx(grid, ind);

    if (bp.is_subsonic_outlet_x() or bp.is_supersonic_outlet_x()) {
        dTxy_dx = 0.0;
        dqx_dx = 0.;
    }
    else if (bp.is_adiabatic_no_slip_wall_x()) {
        auto q_x = [&](
            int i) { return (i == ind) ? 0.0 : der->DX(grid, alias::T, ind); };
        dqx_dx = heat_prefactor * der->DX(grid, q_x, ind);
    }
    auto point = grid.values(ind);
    double e_diss = der->DX(grid, alias::U, ind) * Txx + pf.u(point) * dTxx_dx
        + der->DX(grid, alias::V, ind) * Txy + pf.v(point) * dTxy_dx + dqx_dx;
    return {0., dTxx_dx, dTxy_dx, e_diss};
}

Flux Boundary::dissipation_y(
    const CartesianGrid& grid, const BoundaryPoint& bp) const
{
    int ind = bp.ind;
    double Txy = dissipation_tool->Txy(grid, ind);
    double Tyy = dissipation_tool->Tyy(grid, ind);
    double dTxy_dy = dissipation_tool->dTxy_dy(grid, ind);
    double dTyy_dy = dissipation_tool->dTyy_dy(grid, ind);
    double dqy_dy = dissipation_tool->dqy_dy(grid, ind);

    if (bp.is_subsonic_outlet_y() or bp.is_supersonic_outlet_y()) {
        dTxy_dy = 0.0;
        dqy_dy = 0.;
    }
    else if (bp.is_adiabatic_no_slip_wall_y()) {
        auto q_y = [&](
            int i) { return (i == ind) ? 0.0 : der->DY(grid, alias::T, ind); };
        dqy_dy = heat_prefactor * der->DY(grid, q_y, ind);
    }
    auto point = grid.values(ind);
    double e_diss = der->DY(grid, alias::V, ind) * Tyy + pf.v(point) * dTyy_dy
        + der->DY(grid, alias::U, ind) * Txy + pf.u(point) * dTxy_dy + dqy_dy;
    return {0., dTxy_dy, dTyy_dy, e_diss};
}

Boundary::d_array Boundary::d_from_lodi(const CartesianGrid& grid,
    const LodiArray& lodi_a, const BoundaryPoint& bp) const
{
    double d1, d2, d3, d4;
    double c = pf.sound_speed(grid.values(bp.ind));
    double rho = pf.rho(grid.values(bp.ind));

    d1 = 1 / (c * c) * (lodi_a[1] + 1 / 2. * (lodi_a[0] + lodi_a[3]));
    d2 = 1 / 2. * (lodi_a[0] + lodi_a[3]);
    d3 = 1 / (2 * rho * c) * (lodi_a[3] - lodi_a[0]);
    d4 = lodi_a[2];
    return {{d1, d2, d3, d4}};
}

Flux Boundary::flux_from_lodi(const CartesianGrid& grid,
    const LodiArray& lodi_a, const BoundaryConfiguration& bconf) const
{
    auto d = d_from_lodi(grid, lodi_a, bconf.bp);
    Flux convection{};
    auto point = grid.values(bconf.bp.ind);
    double rho = pf.rho(point);
    double v1 = (pf.*(bconf.v1))(point);
    double v2 = (pf.*(bconf.v2))(point);

    convection.rho = -d[0];
    convection.ru = -v1 * d[0] - rho * d[2];
    convection.rv = -v2 * d[0] - rho * d[3];
    convection.e = -1 / 2. * (v1 * v1 + v2 * v2) * d[0] - d[1] / (gam - 1)
        - rho * v1 * d[2] - rho * v2 * d[3];
    return convection;
}
