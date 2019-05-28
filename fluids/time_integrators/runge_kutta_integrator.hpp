#ifndef RUNGE_KUTTA_INTEGRATOR_HPP
#define RUNGE_KUTTA_INTEGRATOR_HPP

#include "../input_output/writers/writer.hpp"
#include "../utils/point_functions.hpp"
#include "../utils/filters/minimal_filter.hpp"
#include "../utils/filters/minimal_filter_factory.hpp"
#include "../utils/global_vars.hpp"
#include "time_integrator_tool.hpp"
#include <cstdio>
#include <iostream>
#include <memory>
#include <utility>

template <typename Grid, typename Variation>
class RungeKuttaIntegrator {
public:
    RungeKuttaIntegrator(Options& opt_in, Grid& grid_in,
        std::shared_ptr<TimeIntegratorTool> tool_in, PointFunctions& pf_in);
    void run();

private:
    std::shared_ptr<TimeIntegratorTool> tool;
    Grid& grid;
    PointFunctions& pf;
    const double initial_time;
    const double final_time;
    const double cfl;
    Writer writer;
    Grid aux_grid;
    const double print_interval;
    const double reynolds;
    const bool should_filter;
    std::shared_ptr<MinimalFilter> minimal_filter;
    bool found_nan;

    double get_dt(const CartesianGrid& grid_in);
    double get_dt(const ShockGrid& grid_in);
};

template <typename Grid, typename Variation>
RungeKuttaIntegrator<Grid, Variation>::RungeKuttaIntegrator(Options& opt_in,
    Grid& grid_in, std::shared_ptr<TimeIntegratorTool> tool_in,
    PointFunctions& pf_in)
    : tool(std::move(tool_in))
    , grid(grid_in)
    , pf(pf_in)
    , initial_time(opt_in.t_init())
    , final_time(opt_in.t_max())
    , cfl(opt_in.cfl())
    , writer(Writer(opt_in))
    , aux_grid(opt_in)
    , print_interval(opt_in.print_interval())
    , reynolds(opt_in.reynolds())
    , should_filter(opt_in.should_filter())
    , minimal_filter(create_minimal_filter(opt_in.filter_order()))
    , found_nan(false)
{
}

template <typename Grid, typename Variation>
void RungeKuttaIntegrator<Grid, Variation>::run()
{
    double dt;
    double next_print = initial_time + print_interval;
    Variation k1(grid.nPointsTotal);
    Variation k2(grid.nPointsTotal);
    Variation k3(grid.nPointsTotal);
    writer.write(grid, initial_time);
    grid.specific_print();
    for (double t = initial_time; t < final_time; t += dt) {
        dt = get_dt(grid);
        if (found_nan) {
            break;
        }
        if (t + dt >= next_print) {
            dt = next_print - t;
        }
        printf("\r\bt=%e dt=%e", t, dt);
        fflush(nullptr);
        grid.grid_specific_pre_update(dt);
        if (should_filter) {
            minimal_filter->filter_grid(&grid);
        }
        aux_grid.update_values(&grid);
        // Compute k1
        tool->time_derivative(k1, grid, t);
// Compute k2
#ifndef DEBUG
#pragma omp parallel for
#endif
        for (int ind = 0; ind < grid.nPointsTotal; ind++) {
            aux_grid.set_values(
                grid.values(ind) + (dt / 2.) * k1.grid_variation[ind], ind);
        }
        tool->update_values(&aux_grid, t + dt / 2);
        tool->time_derivative(k2, aux_grid, t + dt / 2);
// Compute k3
#ifndef DEBUG
#pragma omp parallel for
#endif
        for (int ind = 0; ind < grid.nPointsTotal; ind++) {
            aux_grid.set_values(grid.values(ind)
                    + (-dt) * k1.grid_variation[ind]
                    + (2 * dt) * k2.grid_variation[ind],
                ind);
        }
        tool->update_values(&aux_grid, t + dt);
        tool->time_derivative(k3, aux_grid, t + dt);
// Update grid
#ifndef DEBUG
#pragma omp parallel for
#endif
        for (int ind = 0; ind < grid.nPointsTotal; ind++) {
            grid.set_values(grid.values(ind)
                    + (1 / 6. * dt) * k1.grid_variation[ind]
                    + (4 / 6. * dt) * k2.grid_variation[ind]
                    + (1 / 6. * dt) * k3.grid_variation[ind],
                ind);
        }
        tool->update_values(&grid, t + dt);
        grid.grid_specific_pos_update(dt);
        if (t + dt >= next_print) {
            writer.write(grid, t + dt);
            grid.specific_print();
            next_print += print_interval;
        }
    }
    writer.write(grid, final_time);
    std::cout << "Exit runge" << std::endl;
}

template <typename Grid, typename Variation>
double RungeKuttaIntegrator<Grid, Variation>::get_dt(
    const CartesianGrid& grid_in)
{
    double dt = 1e10;
    max_mach_number = max_u_plus_c = max_v_plus_c = 0.0;
#ifndef DEBUG
#pragma omp parallel for reduction(min                                         \
                                   : dt),                                      \
    reduction(max                                                              \
              : max_mach_number, max_u_plus_c, max_v_plus_c)
#endif
    for (int ind = 0; ind < grid_in.nPointsTotal; ind++) {
        auto point = grid.values(ind);
        double c = pf.sound_speed(point);
        double u = fabs(pf.u(point));
        double v = fabs(pf.v(point));
        const double dx = grid.dx;
        const double dy = grid.dy;
        double new_dt = 1 / 2.
            * std::min({dx / (u + c), dy / (v + c), dx * dx * reynolds / 2.,
                  dy * dy * reynolds / 2.});
        if (new_dt < dt) {
            dt = new_dt;
        }
        auto new_mach = pf.mach_number(point);
        if (max_mach_number < new_mach) {
            max_mach_number = new_mach;
        }
        auto new_u_plus_c = u + c;
        if (max_u_plus_c < new_u_plus_c) {
            max_u_plus_c = new_u_plus_c;
        }
        auto new_v_plus_c = v + c;
        if (max_v_plus_c < new_v_plus_c) {
            max_v_plus_c = new_v_plus_c;
        }
        if (point.rho() != point.rho() || point.ru() != point.ru()
            || point.rv() != point.rv() || point.e() != point.e()) {
            found_nan = true;
            std::cout << std::endl
                      << "NaN found at ind=" << ind << " x=" << grid.X(ind)
                      << " y=" << grid.Y(ind) << std::endl;
        }
    }
    if (max_mach_number > 1) {
        max_mach_number = 1;
    }
    return cfl * dt;
}

template <typename Grid, typename Variation>
double RungeKuttaIntegrator<Grid, Variation>::get_dt(const ShockGrid& grid_in)
{
    auto dt = get_dt(dynamic_cast<const CartesianGrid&>(grid_in));
    for (auto& sp : grid.shock_points()) {
        double shock_dt;
        if (sp.is_x()) {
            shock_dt = fabs(grid.dx * cos(sp.theta) / sp.w);
        }
        else {
            shock_dt = fabs(grid.dy * sin(sp.theta) / sp.w);
        }
        if (dt > shock_dt) {
            dt = shock_dt * cfl;
        }
    }
    return dt;
}
#endif /* RUNGE_KUTTA_INTEGRATOR_HPP */
