#ifndef EULER_INTEGRATOR_HPP
#define EULER_INTEGRATOR_HPP

#include "../input_output/writers/writer.hpp"
#include "../utils/point_functions.hpp"
#include "time_integrator_tool.hpp"
#include "time_integrator_types.hpp"
#include <algorithm>
#include <iostream>
#include <memory>
#include <utility>

template <typename Grid, typename Variation>
class EulerIntegrator {
public:
    EulerIntegrator(Options& opt_in, Grid& grid_in,
        std::shared_ptr<TimeIntegratorTool> tool_in, PointFunctions& pf_in);
    void run();

private:
    std::shared_ptr<TimeIntegratorTool> tool;
    Grid& grid;
    Variation var;
    PointFunctions& pf;
    const double initial_time;
    const double final_time;
    const double cfl;
    Writer writer;
    const double print_interval;
    const double reynolds;

    double get_dt(const CartesianGrid& grid_in);
};

template <typename Grid, typename Variation>
EulerIntegrator<Grid, Variation>::EulerIntegrator(Options& opt_in,
    Grid& grid_in, std::shared_ptr<TimeIntegratorTool> tool_in,
    PointFunctions& pf_in)
    : tool(std::move(tool_in))
    , grid(grid_in)
    , var(Variation(grid.nPointsTotal))
    , pf(pf_in)
    , initial_time(opt_in.t_init())
    , final_time(opt_in.t_max())
    , cfl(opt_in.cfl())
    , writer(Writer(opt_in))
    , print_interval(opt_in.print_interval())
    , reynolds(opt_in.reynolds())
{
}

template <typename Grid, typename Variation>
void EulerIntegrator<Grid, Variation>::run()
{
    double dt;
    double next_print = initial_time + print_interval;
    Variation k(grid.nPointsTotal);
    writer.write(grid, initial_time);
    for (double t = initial_time; t < final_time; t += dt) {
        dt = get_dt(grid);
        std::cout << std::scientific << "t=" << t << " dt=" << dt << std::endl;
        fflush(nullptr);
        tool->time_derivative(k, grid, t);
        for (int ind = 0; ind < grid.nPointsTotal; ind++) {
            grid.set_values(grid.values(ind) + dt * k.grid_variation[ind], ind);
        }
        tool->update_values(&grid, t);
        if (t + dt >= next_print) {
            writer.write(grid, t);
            next_print += print_interval;
        }
    }
}

template <typename Grid, typename Variation>
double EulerIntegrator<Grid, Variation>::get_dt(const CartesianGrid& grid_in)
{
    double dt = 1e10;
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
    }
    return cfl * dt;
}

#endif /* EULER_INTEGRATOR_HPP */
