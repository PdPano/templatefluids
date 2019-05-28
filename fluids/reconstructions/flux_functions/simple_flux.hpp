#ifndef SIMPLE_FLUX_HPP
#define SIMPLE_FLUX_HPP

#include "flux_interface.hpp"

class SimpleFlux : public FluxInterface {
public:
    SimpleFlux(PointFunctions& pf_in, Options& opt_in);
    Flux fluxX(const Point& p) const override;
    Flux fluxXPositive(const Point& p) const override;
    Flux fluxXNegative(const Point& p) const override;

    Flux fluxY(const Point& p) const override;
    Flux fluxYPositive(const Point& p) const override;
    Flux fluxYNegative(const Point& p) const override;
};

#endif /* SIMPLE_FLUX_HPP */
