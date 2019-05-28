#ifndef DISSIPATION_FACTORY_HPP
#define DISSIPATION_FACTORY_HPP

#include <memory>

class Dissipation;
class Options;
class Derivatives;
struct PointFunctions;

std::shared_ptr<Dissipation> create_dissipation(
    Options& opt, PointFunctions& pf, std::shared_ptr<Derivatives> der);

#endif /* DISSIPATION_FACTORY_HPP */
