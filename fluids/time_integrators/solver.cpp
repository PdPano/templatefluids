#include "solver.hpp"
#include "../grid/cartesian_grid.hpp"
#include "../grid/ghias_grid.hpp"
#include "../grid/ghias_shock_grid.hpp"
#include "../grid/karagiozis_grid.hpp"
#include "../grid/shock_grid.hpp"
#include "../input_output/options.hpp"
#include "../utils/operators_overloads.hpp"
#include "euler_integrator.hpp"
#include "omp.h"
#include "runge_kutta_integrator.hpp"
#include "time_integrator_tool.hpp"
#include "time_integrator_tool_factory.hpp"
#include <memory>

Solver::Solver(Options& opt_in)
    : opt(opt_in)
    , pf(PointFunctions(opt.mach(), opt.gam()))
{
}

void Solver::run()
{
#ifndef DEBUG
    if (opt.omp_threads() > 0) {
        omp_set_num_threads(opt.omp_threads());
    }
#endif
    if (opt.solver_type() == "SIMPLE") {
        setup_and_run<CartesianGrid>();
    }
    else if (opt.solver_type() == "GHIAS") {
        setup_and_run<GhiasGrid>();
    }
    else if (opt.solver_type() == "KARAGIOZIS") {
        setup_and_run<KaragiozisGrid>();
    }
    else if (opt.solver_type() == "SHOCK") {
        setup_and_run<ShockGrid>();
    }
    else if (opt.solver_type() == "GHIAS_SHOCK") {
        setup_and_run<GhiasShockGrid>();
    }
    else {
        std::cout << "Solver " << opt.solver_type() << " not found!!!"
                  << std::endl;
    }
}

template <typename Grid>
void Solver::setup_and_run()
{
    Grid grid(opt);
    auto tool = create_time_integrator_tool(opt, pf, grid);
    if (tool == nullptr) {
        return;
    }
    if (opt.integrator_type() == "EULER") {
        auto integrator
            = EulerIntegrator<Grid, CartesianVariation>(opt, grid, tool, pf);
        integrator.run();
    }
    if (opt.integrator_type() == "RUNGE_KUTTA") {
        auto integrator = RungeKuttaIntegrator<Grid, CartesianVariation>(
            opt, grid, tool, pf);
        integrator.run();
    }
}
