#include "operators_overloads.hpp"
#include "flux_def.hpp"
#include "point_def.hpp"
#include "shock_discontinuity_def.hpp"

Point operator+(const Point& p1, const Point& p2)
{
    return {p1.rho() + p2.rho(), p1.ru() + p2.ru(), p1.rv() + p2.rv(),
        p1.e() + p2.e()};
}

Point operator-(const Point& p1, const Point& p2)
{
    return {p1.rho() - p2.rho(), p1.ru() - p2.ru(), p1.rv() - p2.rv(),
        p1.e() - p2.e()};
}

Point operator/(const Point& p1, const double& d)
{
    return {p1.rho() / d, p1.ru() / d, p1.rv() / d, p1.e() / d};
}

Point operator*(const double& d, const Point& p1)
{
    return {p1.rho() * d, p1.ru() * d, p1.rv() * d, p1.e() * d};
}

Point operator*(const Point& p1, const double& d)
{
    return {p1.rho() * d, p1.ru() * d, p1.rv() * d, p1.e() * d};
}

Flux operator+(const Flux& f1, const Point& p2)
{

    return {f1.rho + p2.rho(), f1.ru + p2.ru(), f1.rv + p2.rv(), f1.e + p2.e()};
}

Point operator+(const Point& p1, const Flux& f2)
{
    return {f2.rho + p1.rho(), f2.ru + p1.ru(), f2.rv + p1.rv(), f2.e + p1.e()};
}

Flux operator*(const Flux& f1, const double& d)
{
    return {f1.rho * d, f1.ru * d, f1.rv * d, f1.e * d};
}

Flux operator*(const double& d, const Flux& f1) { return f1 * d; }

Flux operator-(const Flux& f1) { return {-f1.rho, -f1.ru, -f1.rv, -f1.e}; }

Flux operator/(const Flux& f1, const double& d)
{
    return {f1.rho / d, f1.ru / d, f1.rv / d, f1.e / d};
}

bool operator==(const Point& lhs, const Point& rhs)
{
    return (lhs.rho() == rhs.rho() and lhs.ru() == rhs.ru()
        and lhs.rv() == rhs.rv() and lhs.e() == rhs.e());
}

bool operator==(const ShockDiscontinuity& lhs, const ShockDiscontinuity& rhs)
{
    return (lhs.cleft() == rhs.cleft() and lhs.cright() == rhs.cright()
        and lhs.ind == rhs.ind and lhs.frac == rhs.frac
        and lhs.type == rhs.type);
}
