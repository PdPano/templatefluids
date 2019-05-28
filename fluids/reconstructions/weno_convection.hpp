#ifndef WENO_CONVECTION_HPP
#define WENO_CONVECTION_HPP

#include "abstract_convection.hpp"
#include <memory>
#include <array>

class FluxInterface;
struct Point;

class WenoConvection : public Convection {
public:
    WenoConvection(PointFunctions& pf_in,
        std::shared_ptr<FluxInterface> flux_in,
        std::shared_ptr<Convection> backup_convection_in);
    Flux convection_x(const CartesianGrid& grid, int ind) const override;
    Flux convection_y(const CartesianGrid& grid, int ind) const override;
    void init(const CartesianGrid&) override {}

private:
    using Component = std::array<double, 5>;
    using ComponentsArray = std::array<Component, 4>;
    using Omegas = std::array<double, 3>;
    using OmegasArray = std::array<Omegas, 4>;
    using Betas = std::array<double, 3>;
    using BetasArray = std::array<Betas, 4>;
    using Sigmas = std::array<double, 3>;
    using SigmasArray = std::array<Sigmas, 4>;
    using Approximation = std::array<double, 3>;
    using ApproximationArray = std::array<Approximation, 4>;
    using fluxFunc = Flux (FluxInterface::*)(const Point& p) const;

    std::shared_ptr<FluxInterface> flux;
    std::shared_ptr<Convection> backup_convection;
    Flux flux_at_i_half(
        const CartesianGrid& grid, int ind, fluxFunc func, int shift) const;
    BetasArray compute_betas(const ComponentsArray& vals) const;
    SigmasArray compute_sigmas(const BetasArray& betas) const;
    OmegasArray compute_omegas(const SigmasArray& sigmas) const;
    Omegas remap_omegas(Sigmas sigmas, Omegas omegas);
    ComponentsArray build_components_array(
        const CartesianGrid& grid, int ind, fluxFunc func, int shift) const;
    ApproximationArray compute_approximations(
        const ComponentsArray& vals) const;
};

#endif /* WENO_CONVECTION_HPP */
