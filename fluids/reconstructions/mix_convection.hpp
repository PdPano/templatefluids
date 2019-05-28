#ifndef MIX_CONVECTION_HPP
#define MIX_CONVECTION_HPP

#include "abstract_convection.hpp"
#include <memory>

class FluxInterface;

class MixConvection : public Convection {
public:
    MixConvection(PointFunctions& pf_in, std::shared_ptr<Derivatives> der_in,
        std::shared_ptr<Convection> main_in, std::shared_ptr<Convection> aux_in,
        double mix_param_in);
    Flux convection_x(const CartesianGrid& grid, int ind) const override;
    Flux convection_y(const CartesianGrid& grid, int ind) const override;
    void init(const CartesianGrid&) override {}

private:
    std::shared_ptr<Convection> main;
    std::shared_ptr<Convection> aux;
    const double mix_param;
};

#endif /* MIX_CONVECTION_HPP */
