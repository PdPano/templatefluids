#include "LF_flux.hpp"
#include "../../utils/global_vars.hpp"
#include "../../utils/operators_overloads.hpp"

LaxFriedrichsFlux::LaxFriedrichsFlux(PointFunctions& pf_in, Options& opt_in)
    : FluxInterface(pf_in, opt_in)
{
}

Flux LaxFriedrichsFlux::fluxX(const Point& p) const
{
    return {pf.ru(p), pf.ru2(p) + pf.pressure(p), pf.ruv(p),
        pf.u(p) * (pf.e(p) + pf.pressure(p))};
}

Flux LaxFriedrichsFlux::fluxXPositive(const Point& p) const
{
    return 0.5 * (this->fluxX(p) + max_u_plus_c * p);
}

Flux LaxFriedrichsFlux::fluxXNegative(const Point& p) const
{
    return 0.5 * (this->fluxX(p) + (-1) * max_u_plus_c * p);
}

Flux LaxFriedrichsFlux::fluxY(const Point& p) const
{
    return {pf.rv(p), pf.ruv(p), pf.rv2(p) + pf.pressure(p),
        pf.v(p) * (pf.e(p) + pf.pressure(p))};
}

Flux LaxFriedrichsFlux::fluxYPositive(const Point& p) const
{
    return 0.5 * (this->fluxY(p) + max_v_plus_c * p);
}
Flux LaxFriedrichsFlux::fluxYNegative(const Point& p) const
{
    return 0.5 * (this->fluxY(p) + (-1) * max_v_plus_c * p);
}
