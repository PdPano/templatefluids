#ifndef STEGER_WARMING_FLUX_HPP
#define STEGER_WARMING_FLUX_HPP

#include "flux_interface.hpp"

class StegerWarmingFlux : public FluxInterface {
public:
    StegerWarmingFlux(PointFunctions& pf_in, Options& opt_in);
    Flux fluxX(const Point& p) const override;
    Flux fluxXPositive(const Point& p) const override;
    Flux fluxXNegative(const Point& p) const override;

    Flux fluxY(const Point& p) const override;
    Flux fluxYPositive(const Point& p) const override;
    Flux fluxYNegative(const Point& p) const override;

private:
    template <int k1, int k2, int sign>
    Flux stegerWarmingGenericFlux(const Point& p) const;
    double mixingParameter(double lamb) const;
};

#endif /* STEGER_WARMING_FLUX_HPP */
