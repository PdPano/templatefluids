#ifndef TIME_INTEGRATOR_TOOL_HPP
#define TIME_INTEGRATOR_TOOL_HPP

#include <memory>

class Boundary;
class Dissipation;
class Convection;
class CartesianGrid;
class GhiasShockGrid;
class KaragiozisGrid;
struct CartesianVariation;

class TimeIntegratorTool {
public:
    TimeIntegratorTool(std::shared_ptr<Convection> convection_in,
        std::shared_ptr<Dissipation> dissipation_in,
        std::shared_ptr<Boundary> boundary_in);
    TimeIntegratorTool(std::shared_ptr<Convection> convection_in,
        std::shared_ptr<Dissipation> dissipation_in,
        std::shared_ptr<Boundary> boundary_in,
        std::shared_ptr<Convection> convection_irreg_in,
        std::shared_ptr<Dissipation> dissipation_irreg_in,
        std::shared_ptr<Boundary> boundary_irreg_in);
    void time_derivative(
        CartesianVariation& var, const CartesianGrid& grid, double t);
    void time_derivative(
        CartesianVariation& var, const KaragiozisGrid& grid, double t);
    void time_derivative(
        CartesianVariation& var, const GhiasShockGrid& grid, double t);
    void fix_boundary(CartesianGrid* grid, double t);
    void update_values(CartesianGrid* grid, double t);

private:
    std::shared_ptr<Convection> conv;
    std::shared_ptr<Dissipation> diss;
    std::shared_ptr<Boundary> boundary;

    std::shared_ptr<Convection> conv_irreg;
    std::shared_ptr<Dissipation> diss_irreg;
    std::shared_ptr<Boundary> boundary_irreg;
};

#endif /* TIME_INTEGRATOR_TOOL_HPP */
