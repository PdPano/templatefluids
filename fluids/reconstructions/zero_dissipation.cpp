#include "zero_dissipation.hpp"
#include "../utils/flux_def.hpp"

#include <utility>

ZeroDissipation::ZeroDissipation(PointFunctions& pf_in,
    std::shared_ptr<Derivatives> der_in, double reynolds_in, double prandtl_in)
    : Dissipation(pf_in, std::move(der_in), reynolds_in, prandtl_in)
{
}

Flux ZeroDissipation::dissipation_x(
    const CartesianGrid& /*grid*/, int /*ind*/) const
{
    return {0, 0, 0, 0};
}

Flux ZeroDissipation::dissipation_y(
    const CartesianGrid& /*grid*/, int /*ind*/) const
{
    return {0, 0, 0, 0};
}
