#ifndef LUISA_DETECTOR_FACTORY_HPP
#define LUISA_DETECTOR_FACTORY_HPP

#include <memory>
#include <string>

class Options;
class LuisaDetector;

std::shared_ptr<LuisaDetector> create_luisa_detector(Options& opt, int nPointsI,
    int nPointsJ, std::string detector_type_override = "NONE");

#endif /* LUISA_DETECTOR_FACTORY_HPP */
