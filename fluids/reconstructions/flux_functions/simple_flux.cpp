#include "simple_flux.hpp"
#include "../../utils/operators_overloads.hpp"

SimpleFlux::SimpleFlux(PointFunctions& pf_in, Options& opt_in)
    : FluxInterface(pf_in, opt_in)
{
}

Flux SimpleFlux::fluxX(const Point& p) const
{
    return {pf.ru(p), pf.ru2(p) + pf.pressure(p), pf.ruv(p),
        pf.u(p) * (pf.e(p) + pf.pressure(p))};
}

Flux SimpleFlux::fluxXPositive(const Point& p) const
{
    return 0.5 * (this->fluxX(p));
}

Flux SimpleFlux::fluxXNegative(const Point& p) const
{
    return 0.5 * (this->fluxX(p));
}

Flux SimpleFlux::fluxY(const Point& p) const
{
    return {pf.rv(p), pf.ruv(p), pf.rv2(p) + pf.pressure(p),
        pf.v(p) * (pf.e(p) + pf.pressure(p))};
}

Flux SimpleFlux::fluxYPositive(const Point& p) const
{
    return 0.5 * (this->fluxY(p));
}
Flux SimpleFlux::fluxYNegative(const Point& p) const
{
    return 0.5 * (this->fluxY(p));
}
