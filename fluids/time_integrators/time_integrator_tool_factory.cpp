#include "../boundary/boundary.hpp"
#include "../derivatives/derivatives_factory.hpp"
#include "../input_output/options.hpp"
#include "../reconstructions/abstract_convection.hpp"
#include "../reconstructions/abstract_dissipation.hpp"
#include "../reconstructions/convection_factory.hpp"
#include "../reconstructions/dissipation_factory.hpp"
#include "../utils/operators_overloads.hpp"
#include "time_integrator_tool.hpp"
#include "time_integrator_tool_factory.hpp"
#include <iostream>

std::shared_ptr<TimeIntegratorTool> create_time_integrator_tool(
    Options& opt, PointFunctions& pf, CartesianGrid& grid)
{

    std::shared_ptr<Derivatives> der;
    std::shared_ptr<Convection> conv;
    std::shared_ptr<Dissipation> diss;
    std::shared_ptr<Boundary> boundary;

    der = create_derivative("REGULAR", pf, grid, opt.derivative_order());
    conv = create_convection(opt, pf, der);
    diss = create_dissipation(opt, pf, der);
    boundary
        = std::make_shared<Boundary>(pf, der, opt.reynolds(), opt.prandtl());

    if (der == nullptr or conv == nullptr or diss == nullptr) {
        std::cerr << "Could not create tool!!!" << std::endl;
        return nullptr;
    }

    if (opt.solver_type() == "KARAGIOZIS" || opt.solver_type() == "SHOCK") {
        auto der_irreg = create_derivative("KARAGIOZIS", pf, grid);
        auto conv_irreg = create_convection(opt, pf, der_irreg);
        auto diss_irreg = create_dissipation(opt, pf, der_irreg);
        auto boundary_irreg = std::make_shared<Boundary>(
            pf, der_irreg, opt.reynolds(), opt.prandtl());
        if (der_irreg == nullptr or conv_irreg == nullptr
            or diss_irreg == nullptr) {
            std::cerr << "Could not create tool!!! (Irregular part)"
                      << std::endl;
            return nullptr;
        }
        return std::make_shared<TimeIntegratorTool>(
            conv, diss, boundary, conv_irreg, diss_irreg, boundary_irreg);
    }

    if (opt.solver_type() == "GHIAS_SHOCK") {
        auto conv_irreg = create_convection(opt, pf, der, "WENO_CONVECTION");
        auto diss_irreg = create_dissipation(opt, pf, der);
        return std::make_shared<TimeIntegratorTool>(
            conv, diss, boundary, conv_irreg, diss, boundary);
    }

    return std::make_shared<TimeIntegratorTool>(conv, diss, boundary);
}
