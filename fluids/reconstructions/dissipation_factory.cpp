#include "dissipation_factory.hpp"
#include "../input_output/options.hpp"
#include "../utils/point_functions.hpp"
#include "abstract_dissipation.hpp"
#include "simple_dissipation.hpp"
#include "zero_dissipation.hpp"

std::shared_ptr<Dissipation> create_dissipation(
    Options& opt, PointFunctions& pf, std::shared_ptr<Derivatives> der)
{
    if (opt.dissipation() == "SIMPLE") {
        return std::make_shared<SimpleDissipation>(
            pf, der, opt.reynolds(), opt.prandtl());
    }
    if (opt.dissipation() == "ZERO") {
        return std::make_shared<ZeroDissipation>(
            pf, der, opt.reynolds(), opt.prandtl());
    }
    return (nullptr);
}
