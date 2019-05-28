#include "../../derivatives/derivatives.hpp"
#include "../../input_output/options.hpp"
#include "../../utils/point_functions.hpp"
#include "LF_flux.hpp"
#include "flux_interface.hpp"
#include "simple_flux.hpp"
#include "steger_warming_flux.hpp"
#include <iostream>

std::shared_ptr<FluxInterface> create_flux(
    Options& opt, PointFunctions& pf, std::string flux_type_override)
{
    if (flux_type_override == "NONE") { // NO_OVERRIDE
        flux_type_override = opt.flux();
    }
    if (flux_type_override == "SIMPLE") {
        return std::make_shared<SimpleFlux>(pf, opt);
    }
    if (flux_type_override == "STEGER_WARMING") {
        return std::make_shared<StegerWarmingFlux>(pf, opt);
    }
    if (flux_type_override == "LF_FLUX") {
        return std::make_shared<LaxFriedrichsFlux>(pf, opt);
    }
    std::cerr << "Flux type " << flux_type_override
              << " not found! Using SIMPLE instead" << std::endl;
    return std::make_shared<SimpleFlux>(pf, opt);
}
