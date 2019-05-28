#ifndef CONVECTION_FACTORY_HPP
#define CONVECTION_FACTORY_HPP

#include <memory>

class Convection;
class Derivatives;
class Options;
struct PointFunctions;

std::shared_ptr<Convection> create_convection(Options& opt, PointFunctions& pf,
    std::shared_ptr<Derivatives> der, std::string overwrite_conv = "NONE");

#endif /* CONVECTION_FACTORY_HPP */
