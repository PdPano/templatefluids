#include "shock_discontinuity_handler.hpp"
#include "../input_output/options.hpp"
#include "../utils/debug_def.hpp"
#include "point_def.hpp"
#include "point_functions.hpp"
#include <iostream>

ShockHandler::ShockHandler(Options& opt)
    : pf(PointFunctions(opt.mach(), opt.gam()))
    , gam(opt.gam())
    , delta((gam - 1.) / 2)
    , angular_coeff(3.
          / (sqrt((16. * gam - delta) * (1 + 16 * delta)) / (4 * (1. + delta))
                + 15. * delta / (4 * (1. + delta)) - 1))
    , cfl(opt.cfl())
{
}

ShockDiscontinuity ShockHandler::create_shock(
    Point left_in, Point right_in, char type_in, double dx_dt)
{
    ShockDiscontinuity shock(left_in, right_in, type_in);
    shock.low_pressure_side
        = (pf.pressure(shock.left()) < pf.pressure(shock.right())) ? 'l' : 'r';
    fix_theta(shock);
    if (is_candidate(shock, dx_dt)) {
        prepare_values(shock);
    }
    return shock;
}

void ShockHandler::prepare_values(ShockDiscontinuity& shock)
{
    ShockHandler::Velocity low_vel
        = to_normal(vel_from_point(low_point(shock)), shock.theta);
    ShockHandler::Velocity high_vel
        = to_normal(vel_from_point(high_point(shock)), shock.theta);
    shock.sigma = (high_c(shock) + delta * fabs(u(low_vel) - u(high_vel)))
        / low_c(shock);

    double rel_mach = calculate_rel_mach(shock.sigma);
    double w = calculate_w(rel_mach, u(low_vel), low_c(shock));
    shock.w = w;

    apply_rh_from_w_mach(shock, low_vel, w, rel_mach);
}

void ShockHandler::compute_wall_interaction(
    ShockDiscontinuity& shock, double wall_theta)
{
    // Collision with wall change which side has low pressure
    shock.low_pressure_side = (shock.low_pressure_side == 'l') ? 'r' : 'l';

    // Shock angle is that of the wall
    shock.theta = wall_theta;
    fix_theta(shock);

    ShockHandler::Velocity low_vel
        = to_normal(vel_from_point(low_point(shock)), shock.theta);
    double rel_mach = get_mach_at_wall_col(u(low_vel), low_c(shock));
    double w = calculate_w(rel_mach, u(low_vel), low_c(shock));
    shock.w = w;
    apply_rh_from_w_mach(shock, low_vel, w, rel_mach);
    shock.sigma = sigma_from_mach(rel_mach);
}

double ShockHandler::get_mach_at_wall_col(double low_u, double low_c)
{
    return low_u * (1 + delta) / (2 * low_c)
        + sqrt(pow(low_u * (1 + delta) / (2 * low_c), 2) + 1);
}

void ShockHandler::apply_rh_from_w_mach(ShockDiscontinuity& shock,
    ShockHandler::Velocity& low_vel, double w, double rel_mach)
{
    /*Rankine-Hugoniot jump conditions*/
    double rel_mach2 = rel_mach * rel_mach;
    double u_high = u(low_vel)
        + low_c(shock) * (1.0 - rel_mach2) / ((1 + delta) * rel_mach);
    double v_high = v(low_vel);
    double rho_low = pf.rho(low_point(shock));
    double rho_high = rho_low * (u(low_vel) - w) / (u_high - w);
    double c_high = low_c(shock)
        * sqrt((gam * rel_mach2 - delta) * (1 + delta * rel_mach2))
        / ((1 + delta) * rel_mach);
    double p_high = c_high * c_high * rho_high / gam;
    double E_high = p_high / (gam - 1)
        + 1 / 2. * rho_high * (u_high * u_high + v_high * v_high);
    auto high_vel = from_normal(std::make_pair(u_high, v_high), shock.theta);
    high_point(shock) = Point(
        rho_high, rho_high * u(high_vel), rho_high * v(high_vel), E_high);
}

void ShockHandler::fix_theta(ShockDiscontinuity& shock)
{
    double theta = shock.theta;
    if (shock.is_x()) {
        if (-M_PI / 2 < theta and theta <= M_PI / 2) {
            if (shock.low_pressure_side == 'l') {
                shock.theta = theta;
            }
            else {
                shock.theta = theta + M_PI;
            }
        }
        else {
            if (shock.low_pressure_side == 'l') {
                shock.theta = theta + M_PI;
            }
            else {
                shock.theta = theta;
            }
        }
    }
    if (shock.is_y()) {
        if (0 < theta and theta <= M_PI) {
            if (shock.low_pressure_side == 'l') {
                shock.theta = theta;
            }
            else {
                shock.theta = theta - M_PI;
            }
        }
        else {
            if (shock.low_pressure_side == 'l') {
                shock.theta = theta + M_PI;
            }
            else {
                shock.theta = theta;
            }
        }
    }
    while (shock.theta > M_PI) {
        shock.theta -= 2 * M_PI;
    }
    while (shock.theta < -M_PI) {
        shock.theta += 2 * M_PI;
    }
}

double ShockHandler::sigma_from_mach(double mach)
{
    return (sqrt((gam * mach * mach - delta) * (1 + delta * mach * mach))
               + delta * (mach * mach - 1))
        / (mach * (1 + delta));
}

Point& ShockHandler::low_point(ShockDiscontinuity& shock)
{
    if (shock.low_pressure_side == 'l') {
        return shock.left();
    }
    return shock.right();
}

Point& ShockHandler::high_point(ShockDiscontinuity& shock)
{
    if (shock.low_pressure_side == 'r') {
        return shock.left();
    }
    return shock.right();
}

double ShockHandler::estimate_rel_mach(double sigma)
{
    return 1. + angular_coeff * (sigma - 1);
}

double ShockHandler::low_c(ShockDiscontinuity& shock)
{
    return pf.sound_speed(low_point(shock));
}

double ShockHandler::high_c(ShockDiscontinuity& shock)
{
    return pf.sound_speed(high_point(shock));
}

ShockHandler::Velocity ShockHandler::to_normal(
    ShockHandler::Velocity vel, double theta)
{
    return std::make_pair(cos(theta) * u(vel) + sin(theta) * v(vel),
        -sin(theta) * u(vel) + cos(theta) * v(vel));
}

ShockHandler::Velocity ShockHandler::from_normal(
    ShockHandler::Velocity vel, double theta)
{
    return to_normal(vel, -theta);
}

double ShockHandler::u(ShockHandler::Velocity vel) { return vel.first; }

double ShockHandler::v(ShockHandler::Velocity vel) { return vel.second; }

ShockHandler::Velocity ShockHandler::vel_from_point(Point& p)
{
    return std::make_pair(pf.u(p), pf.v(p));
}

double ShockHandler::calculate_w(double rel_mach, double low_u, double low_c)
{
    return low_u - rel_mach * low_c;
}

double ShockHandler::calculate_rel_mach(double sigma)
{
    double new_mach = estimate_rel_mach(sigma);
    double estimated_sigma = sigma_from_mach(new_mach);
    double delta_sigma = sigma - estimated_sigma;
    while (fabs(delta_sigma) > 1e-10) {
        new_mach += angular_coeff * delta_sigma;
        estimated_sigma = sigma_from_mach(new_mach);
        delta_sigma = sigma - estimated_sigma;
    }
    return new_mach;
}

bool ShockHandler::is_candidate(ShockDiscontinuity& shock, double dx_dt)
{
    double u_left, u_right;
    if (shock.is_x()) {
        u_left = pf.u(shock.left());
        u_right = pf.u(shock.right());
    }
    else {
        u_left = pf.v(shock.left());
        u_right = pf.v(shock.right());
    }
    double c_left = pf.sound_speed(shock.left());
    double c_right = pf.sound_speed(shock.right());
    double delta_lamb;

    if (shock.low_pressure_side == 'l') {
        delta_lamb = ((u_right - c_right) - (u_left - c_left));
    }
    else {
        delta_lamb = ((u_right + c_right) - (u_left + c_left));
    }

    /*delta_lamb is negative when it is a candidate*/
    if (-delta_lamb / dx_dt > 0.22 * cfl) {
        return true;
    }
    shock.sigma = 1.0;
    return false;
}
