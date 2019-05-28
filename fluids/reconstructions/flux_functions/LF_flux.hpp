#ifndef LF_FLUX_HPP
#define LF_FLUX_HPP

#include "flux_interface.hpp"

class LaxFriedrichsFlux : public FluxInterface {
public:
    LaxFriedrichsFlux(PointFunctions& pf_in, Options& opt_in);
    Flux fluxX(const Point& p) const override;
    Flux fluxXPositive(const Point& p) const override;
    Flux fluxXNegative(const Point& p) const override;

    Flux fluxY(const Point& p) const override;
    Flux fluxYPositive(const Point& p) const override;
    Flux fluxYNegative(const Point& p) const override;
};

#endif /* LF_FLUX_HPP */
