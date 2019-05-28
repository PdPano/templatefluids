#ifndef SIMPLE_DISSIPATION_HPP
#define SIMPLE_DISSIPATION_HPP
#include "abstract_dissipation.hpp"

class SimpleDissipation : public Dissipation {

public:
    SimpleDissipation(PointFunctions& pf_in,
        std::shared_ptr<Derivatives> der_in, double reynolds_in,
        double prandtl_in);
    Flux dissipation_x(const CartesianGrid& grid, int ind) const override;
    Flux dissipation_y(const CartesianGrid& grid, int ind) const override;
};

#endif /* SIMPLE_DISSIPATION_HPP */
