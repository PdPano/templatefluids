#include "luisa_detector_factory.hpp"
#include "../../input_output/options.hpp"
#include "luisa_detector_23.hpp"
#include "luisa_detector_345.hpp"
#include "luisa_shock_detector.hpp"
#include <iostream>

std::shared_ptr<LuisaDetector> create_luisa_detector(Options& opt, int nPointsI,
    int nPointsJ, std::string detector_type_override)
{
    if (detector_type_override == "NONE") { // NO_OVERRIDE
        detector_type_override = opt.luisa_detector();
    }
    if (detector_type_override == "TYPE_23") {
        return std::make_shared<LuisaDetector23>(
            opt.detector_sensitivity(), nPointsI, nPointsJ);
    }
    if (detector_type_override == "TYPE_345") {
        return std::make_shared<LuisaDetector345>(
            opt.detector_sensitivity(), nPointsI, nPointsJ);
    }
    std::cerr << "Detector type " << detector_type_override
              << " not found! Using TYPE_23 instead" << std::endl;
    return std::make_shared<LuisaDetector23>(
        opt.detector_sensitivity(), nPointsI, nPointsJ);
}
