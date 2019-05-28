#include "time_integrator_tool.hpp"
#include "../boundary/boundary.hpp"
#include "../grid/cartesian_grid.hpp"
#include "../grid/ghias_shock_grid.hpp"
#include "../grid/karagiozis_grid.hpp"
#include "../reconstructions/abstract_convection.hpp"
#include "../reconstructions/abstract_dissipation.hpp"
#include "time_integrator_types.hpp"
#include <omp.h>

#include <utility>

TimeIntegratorTool::TimeIntegratorTool(
    std::shared_ptr<Convection> convection_in,
    std::shared_ptr<Dissipation> dissipation_in,
    std::shared_ptr<Boundary> boundary_in)
    : TimeIntegratorTool(std::move(convection_in), std::move(dissipation_in),
          std::move(boundary_in), nullptr, nullptr, nullptr)
{
}

TimeIntegratorTool::TimeIntegratorTool(
    std::shared_ptr<Convection> convection_in,
    std::shared_ptr<Dissipation> dissipation_in,
    std::shared_ptr<Boundary> boundary_in,
    std::shared_ptr<Convection> convection_irreg_in,
    std::shared_ptr<Dissipation> dissipation_irreg_in,
    std::shared_ptr<Boundary> boundary_irreg_in)
    : conv(std::move(convection_in))
    , diss(std::move(dissipation_in))
    , boundary(std::move(boundary_in))
    , conv_irreg(std::move(convection_irreg_in))
    , diss_irreg(std::move(dissipation_irreg_in))
    , boundary_irreg(std::move(boundary_irreg_in))
{
}

void TimeIntegratorTool::time_derivative(
    CartesianVariation& var, const CartesianGrid& grid, double /*t*/)
{
    int nPointsTotal = grid.nPointsTotal;
    conv->init(grid);
#ifndef DEBUG
#pragma omp parallel for
#endif
    for (int ind = 0; ind < nPointsTotal; ind++) { // NOLINT
        var.grid_variation[ind] = conv->convection_x(grid, ind)
            /*+ conv->convection_y(grid, ind)*/
            + diss->dissipation_x(grid, ind)
            /*+ diss->dissipation_y(grid, ind)*/;
    }
    for (auto& bp : grid.boundary()) { // NOLINT
        int ind = bp.ind;
        var.grid_variation[ind] = {0, 0, 0, 0};
        /* if (bp.x_boundary) {
             var.grid_variation[ind] = (boundary->convection_x(grid, bp, t)
                 + boundary->dissipation_x(grid, bp));
         }
         else {
             var.grid_variation[ind] = (conv->convection_x(grid, ind)
                 + diss->dissipation_x(grid, ind));
         }
         if (bp.y_boundary) {
            var.grid_variation[ind] +=  (boundary->convection_y(grid, bp, t)
                  + boundary->dissipation_y(grid, bp));

         }
         else {
             var.grid_variation[ind] += (conv->convection_y(grid, ind)
                  + diss->dissipation_y(grid, ind));
         }*/
    }
}

void TimeIntegratorTool::time_derivative(
    CartesianVariation& var, const KaragiozisGrid& grid, double /*t*/)
{

    int nPointsTotal = grid.nPointsTotal;
#ifndef DEBUG
#pragma omp parallel for
#endif
    for (int ind = 0; ind < nPointsTotal; ind++) { // NOLINT
        var.grid_variation[ind] = conv->convection_x(grid, ind)
            /*+ conv->convection_y(grid, ind)*/
            + diss->dissipation_x(grid, ind)
            /*+ diss->dissipation_y(grid, ind)*/;
    }
    for (const auto& ind : grid.to_revisit()) {
        var.grid_variation[ind] = conv_irreg->convection_x(grid, ind)
            /*+ conv_irreg->convection_y(grid, ind)*/
            + diss_irreg->dissipation_x(grid, ind)
            /*+ diss_irreg->dissipation_y(grid, ind)*/;
    }
    for (auto& bp : grid.boundary()) { // NOLINT
        int ind = bp.ind;
        var.grid_variation[ind] = {0, 0, 0, 0};
        /*if (bp.x_boundary) {
            var.grid_variation[ind] = (boundary_irreg->convection_x(grid, bp, t)
                + boundary_irreg->dissipation_x(grid, bp));
        }
        else {
            var.grid_variation[ind] = (conv_irreg->convection_x(grid, ind)
                + diss_irreg->dissipation_x(grid, ind));
        }
        if (bp.y_boundary) {
            var.grid_variation[ind]
                += (boundary_irreg->convection_y(grid, bp, t)
                     + boundary_irreg->dissipation_y(grid, bp));
        }
        else {
            var.grid_variation[ind]
                += (conv_irreg->convection_y(grid, ind)
+ diss_irreg->dissipation_y(grid, ind));
        }*/
    }
}

void TimeIntegratorTool::time_derivative(
    CartesianVariation& var, const GhiasShockGrid& grid, double /*t*/)
{

    int nPointsTotal = grid.nPointsTotal;
#ifndef DEBUG
#pragma omp parallel for
#endif
    for (int ind = 0; ind < nPointsTotal; ind++) { // NOLINT
        var.grid_variation[ind] = conv->convection_x(grid, ind)
            /*+ conv->convection_y(grid, ind)*/
            + diss->dissipation_x(grid, ind)
            /*+ diss->dissipation_y(grid, ind)*/;
    }

    auto& to_revisit = grid.to_revisit();
    for (auto& ind : to_revisit) {
        var.grid_variation[ind] = conv_irreg->convection_x(grid, ind)
            /*+ conv_irreg->convection_y(grid, ind)*/
            + diss_irreg->dissipation_x(grid, ind)
            /*+ diss_irreg->dissipation_y(grid, ind)*/;
    }
    for (auto& bp : grid.boundary()) { // NOLINT
        int ind = bp.ind;
        var.grid_variation[ind] = {0., 0., 0., 0.};
        /*
        if (bp.x_boundary) {
           var.grid_variation[ind] = (boundary->convection_x(grid, bp, t)
                + boundary->dissipation_x(grid, bp));
        }
        else {
            var.grid_variation[ind] = conv_irreg->convection_x(grid, ind)
                + diss_irreg->dissipation_x(grid, ind);
        }

        if (bp.y_boundary) {
            var.grid_variation[ind]
                +=  (boundary->convection_y(grid, bp, t)
+ boundary->dissipation_y(grid, bp));
        }
        else {
            var.grid_variation[ind]
                += conv_irreg->convection_y(grid, ind)
+ diss_irreg->dissipation_y(grid, ind);
        }*/
    }
}

void TimeIntegratorTool::fix_boundary(CartesianGrid* /*grid*/, double /*t*/)
{
    /*for (auto& bp : grid->boundary()) {
        boundary->fix_boundary(grid, bp, t);
    }*/
}

void TimeIntegratorTool::update_values(CartesianGrid* grid, double t)
{
    fix_boundary(grid, t);
    grid->grid_specific_update();
}
