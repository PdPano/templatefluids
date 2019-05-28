#ifndef SPLIT_CONVECTION_CACHED_HPP
#define SPLIT_CONVECTION_CACHED_HPP

#include "abstract_convection.hpp"
#include <memory>
#include <vector>

class FluxInterface;

class SplitConvectionCached : public Convection {
public:
    SplitConvectionCached(PointFunctions& pf_in,
        std::shared_ptr<Derivatives> der_in,
        std::shared_ptr<FluxInterface> flux_in);
    Flux convection_x(const CartesianGrid& grid, int ind) const override;
    Flux convection_y(const CartesianGrid& grid, int ind) const override;
    void init(const CartesianGrid& grid) override;

private:
    std::shared_ptr<FluxInterface> flux;
    std::vector<Flux> fluxXPos;
    std::vector<Flux> fluxXNeg;
    std::vector<Flux> fluxYPos;
    std::vector<Flux> fluxYNeg;
};

#endif /* SPLIT_CONVECTION_CACHED_HPP */
