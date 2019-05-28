#ifndef DISSIPATION_TOOL_HPP
#define DISSIPATION_TOOL_HPP

#include <memory>
class CartesianGrid;
class Derivatives;
struct PointFunctions;

class DissipationTool {
public:
    DissipationTool(const PointFunctions& pf_in,
        std::shared_ptr<Derivatives> der_in, double reynolds_in,
        double prandtl_in);
    double Txx(const CartesianGrid& grid, int ind) const;
    double Txy(const CartesianGrid& grid, int ind) const;
    double Tyy(const CartesianGrid& grid, int ind) const;

    double dTxx_dx(const CartesianGrid& grid, int ind) const;
    double dTxy_dx(const CartesianGrid& grid, int ind) const;
    double dTxy_dy(const CartesianGrid& grid, int ind) const;
    double dTyy_dy(const CartesianGrid& grid, int ind) const;

    double dqx_dx(const CartesianGrid& grid, int ind) const;
    double dqy_dy(const CartesianGrid& grid, int ind) const;

private:
    std::shared_ptr<Derivatives> der;
    const double gam;
    const double mach;
    const double reynolds;
    const double prandtl;
    const double heat_prefactor;
};

#endif /* DISSIPATION_TOOL_HPP */
