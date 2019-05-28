#ifndef ABSTRACT_DISSIPATION_HPP
#define ABSTRACT_DISSIPATION_HPP

#include "dissipation_tool.hpp"

class CartesianGrid;
class Derivatives;
struct PointFunctions;
struct Flux;

class Dissipation {
public:
    const PointFunctions& pf;
    std::shared_ptr<Derivatives> der;
    const DissipationTool dissipation_tool;

    virtual Flux dissipation_x(const CartesianGrid& grid, int ind) const = 0;
    virtual Flux dissipation_y(const CartesianGrid& grid, int ind) const = 0;

protected:
    Dissipation(PointFunctions& pf_in, std::shared_ptr<Derivatives> der_in,
        double reynolds_in, double prandtl_in)
        : pf(pf_in)
        , der(der_in)
        , dissipation_tool(pf_in, der_in, reynolds_in, prandtl_in)
    {
    }
};

#endif /* ABSTRACT_DISSIPATION_HPP */
