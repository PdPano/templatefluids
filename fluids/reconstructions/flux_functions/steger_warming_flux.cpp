#include "steger_warming_flux.hpp"
#include <cmath>

StegerWarmingFlux::StegerWarmingFlux(PointFunctions& pf_in, Options& opt_in)
    : FluxInterface(pf_in, opt_in)
{
}

Flux StegerWarmingFlux::fluxX(const Point& p) const
{
    return {pf.ru(p), pf.ru2(p) + pf.pressure(p), pf.ruv(p),
        pf.u(p) * (pf.e(p) + pf.pressure(p))};
}

Flux StegerWarmingFlux::fluxXPositive(const Point& p) const
{
    return this->stegerWarmingGenericFlux<1, 0, 1>(p);
}

Flux StegerWarmingFlux::fluxXNegative(const Point& p) const
{
    return this->stegerWarmingGenericFlux<1, 0, -1>(p);
}

Flux StegerWarmingFlux::fluxY(const Point& p) const
{
    return {pf.rv(p), pf.ruv(p), pf.rv2(p) + pf.pressure(p),
        pf.v(p) * (pf.e(p) + pf.pressure(p))};
}

Flux StegerWarmingFlux::fluxYPositive(const Point& p) const
{
    return this->stegerWarmingGenericFlux<0, 1, 1>(p);
}
Flux StegerWarmingFlux::fluxYNegative(const Point& p) const
{
    return this->stegerWarmingGenericFlux<0, 1, -1>(p);
}

template <int k1, int k2, int sign>
Flux StegerWarmingFlux::stegerWarmingGenericFlux(const Point& p) const
{
    double lamb1, lamb2, lamb3, lamb4;
    double rho, u, v, c;

    c = pf.sound_speed(p);
    u = pf.u(p);
    v = pf.v(p);
    rho = pf.rho(p);

    lamb1 = lamb2 = k1 * u + k2 * v;
    lamb3 = lamb1 + c;
    lamb4 = lamb1 - c;
    lamb1 = (lamb1 + sign * fabs(lamb1)) / 2.0 + sign * mixingParameter(lamb1);
    lamb2 = (lamb2 + sign * fabs(lamb2)) / 2.0 + sign * mixingParameter(lamb2);
    lamb3 = (lamb3 + sign * fabs(lamb3)) / 2.0 + sign * mixingParameter(lamb3);
    lamb4 = (lamb4 + sign * fabs(lamb4)) / 2.0 + sign * mixingParameter(lamb4);

    const double ret_rho
        = (rho / (2 * gam)) * (2 * (gam - 1) * lamb1 + lamb3 + lamb4);
    const double ret_ru
        = (rho / (2 * gam)) * (2 * (gam - 1) * lamb1 * u + lamb3 * (u + c * k1)
                                  + lamb4 * (u - c * k1));
    const double ret_rv
        = (rho / (2 * gam)) * (2 * (gam - 1) * lamb1 * v + lamb3 * (v + c * k2)
                                  + lamb4 * (v - c * k2));
    const double ret_e
        = (rho / (2 * gam))
        * ((gam - 1) * lamb1 * (u * u + v * v)
              + 0.5 * lamb3 * (pow(u + c * k1, 2) + pow(v + c * k2, 2))
              + 0.5 * lamb4 * (pow(u - c * k1, 2) + pow(v - c * k2, 2))
              + (3 - gam) * (lamb3 + lamb4) * c * c / (2 * (gam - 1)));

    return {ret_rho, ret_ru, ret_rv, ret_e};
}

double StegerWarmingFlux::mixingParameter(double lamb) const
{
    return 0 * (1e-1) * exp(-1 * fabs(lamb));
}
