#ifndef SIMPLE_CONVECTION_HPP
#define SIMPLE_CONVECTION_HPP

#include "abstract_convection.hpp"

class SimpleConvection : public Convection {
public:
    SimpleConvection(
        PointFunctions& pf_in, std::shared_ptr<Derivatives> der_in);
    Flux convection_x(const CartesianGrid& grid, int ind) const override;
    Flux convection_y(const CartesianGrid& grid, int ind) const override;
    void init(const CartesianGrid&) override {}
};

#endif /* SIMPLE_CONVECTION_HPP */
