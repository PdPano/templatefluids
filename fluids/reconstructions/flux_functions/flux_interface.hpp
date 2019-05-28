#ifndef FLUX_INTERFACE_HPP
#define FLUX_INTERFACE_HPP

#include "../../input_output/options.hpp"
#include "../../utils/flux_def.hpp"
#include "../../utils/point_def.hpp"
#include "../../utils/point_functions.hpp"
class FluxInterface {
public:
    const PointFunctions& pf;
    const double gam, mach;
    virtual Flux fluxX(const Point& p) const = 0;
    virtual Flux fluxXPositive(const Point& p) const = 0;
    virtual Flux fluxXNegative(const Point& p) const = 0;

    virtual Flux fluxY(const Point& p) const = 0;
    virtual Flux fluxYPositive(const Point& p) const = 0;
    virtual Flux fluxYNegative(const Point& p) const = 0;

protected:
    FluxInterface(PointFunctions& pf_in, Options& opt_in)
        : pf(pf_in)
        , gam(opt_in.gam())
        , mach(opt_in.mach())
    {
    }
};

#endif /* FLUX_INTERFACE_HPP */
