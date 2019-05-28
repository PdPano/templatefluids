#ifndef SKEW_SYMMETRIC_HPP
#define SKEW_SYMMETRIC_HPP

#include "abstract_convection.hpp"

class SkewSymmetric : public Convection {
public:
    SkewSymmetric(PointFunctions& pf_in, std::shared_ptr<Derivatives> der_in);
    Flux convection_x(const CartesianGrid& grid, int ind) const override;
    Flux convection_y(const CartesianGrid& grid, int ind) const override;
    void init(const CartesianGrid&) override {}
};

#endif /* SKEW_SYMMETRIC_HPP */
