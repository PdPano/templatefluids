#ifndef ZERO_DISSIPATION_HPP
#define ZERO_DISSIPATION_HPP

#include "abstract_dissipation.hpp"
#include <memory>

class ZeroDissipation : public Dissipation {
public:
    ZeroDissipation(PointFunctions& pf_in, std::shared_ptr<Derivatives> der_in,
        double reynolds_in, double prandtl_in);
    Flux dissipation_x(const CartesianGrid& grid, int ind) const override;
    Flux dissipation_y(const CartesianGrid& grid, int ind) const override;
};

#endif /* ZERO_DISSIPATION_HPP */
