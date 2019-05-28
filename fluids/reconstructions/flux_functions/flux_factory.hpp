
#ifndef FLUX_FACTORY_HPP
#define FLUX_FACTORY_HPP
#include "flux_factory.hpp"
#include "../../input_output/options.hpp"
#include "../../utils/point_functions.hpp"
#include "flux_interface.hpp"
#include "simple_flux.hpp"
#include "steger_warming_flux.hpp"
#include <memory>
#include <string>

std::shared_ptr<FluxInterface> create_flux(
    Options& opt, PointFunctions& pf, std::string flux_type_override = "NONE");

#endif /* FLUX_FACTORY_HPP */
