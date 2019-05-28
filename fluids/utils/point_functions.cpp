#include "point_functions.hpp"
#include "./point_def.hpp"
#include <cmath>

PointFunctions::PointFunctions(double mach_in, double gam_in)
    : mach(mach_in)
    , gam(gam_in)
{
}

double PointFunctions::rho(const Point& p) const { return p.rho(); }

double PointFunctions::ru(const Point& p) const { return p.ru(); }

double PointFunctions::rv(const Point& p) const { return p.rv(); }

double PointFunctions::e(const Point& p) const { return p.e(); }

double PointFunctions::pressure(const Point& p) const
{
    return (gam - 1) * rho_internal_energy(p);
}

double PointFunctions::sound_speed(const Point& p) const
{
    return sqrt(gam * pressure(p) / rho(p));
}

double PointFunctions::temperature(const Point& p) const
{
    return gam * mach * mach * pressure(p) / rho(p);
}

double PointFunctions::ruv(const Point& p) const { return ru(p) * v(p); }

double PointFunctions::u(const Point& p) const { return ru(p) / rho(p); }

double PointFunctions::v(const Point& p) const { return rv(p) / rho(p); }

double PointFunctions::ru2(const Point& p) const
{
    const double ru_v = ru(p);
    return ru_v * ru_v / rho(p);
}

double PointFunctions::rv2(const Point& p) const
{
    const double rv_v = rv(p);
    return rv_v * rv_v / rho(p);
}

double PointFunctions::internal_energy(const Point& p) const
{
    return rho_internal_energy(p) / rho(p);
}

double PointFunctions::rho_internal_energy(const Point& p) const
{
    return e(p) - (ru2(p) + rv2(p)) / 2.0;
}

double PointFunctions::e_from_T(const Point& p, const double& T) const
{
    return T * rho(p) / (gam * (gam - 1) * mach * mach)
        + 1 / 2. * (ru2(p) + rv2(p));
}

double PointFunctions::entropy(const Point& p) const
{
    return 1 / (gam * (gam - 1)) * log(pressure(p) / pow(rho(p), gam));
}

double PointFunctions::mach_number(const Point& p) const
{
    return sqrt(u(p) * u(p) + v(p) * v(p)) / sound_speed(p);
}
