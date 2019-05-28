#ifndef SIMPLE_FLUX_CONVECTION_HPP
#define SIMPLE_FLUX_CONVECTION_HPP

#include "abstract_convection.hpp"
#include <memory>

class FluxInterface;

class SimpleFluxConvection : public Convection {
public:
    SimpleFluxConvection(PointFunctions& pf_in,
        std::shared_ptr<Derivatives> der_in,
        std::shared_ptr<FluxInterface> flux_in);
    Flux convection_x(const CartesianGrid& grid, int ind) const override;
    Flux convection_y(const CartesianGrid& grid, int ind) const override;
    void init(const CartesianGrid&) override {}

private:
    std::shared_ptr<FluxInterface> flux;
};

#endif /* SIMPLE_FLUX_CONVECTION_HPP */
