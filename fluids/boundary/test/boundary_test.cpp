#include "../boundary.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <memory>
#include <sstream>

#include "../../derivatives/regular_derivatives.hpp"
#include "../../grid/test/cartesian_grid_test_interface.hpp"
#include "../../input_output/options.hpp"
#include "../../utils/point_functions.hpp"
#include "sample_inputs.inc"

class BoundaryTest : public testing::Test {
protected:
    virtual void SetUp()
    {

        std::istringstream initial_conditions(initial_conditions_sample);
        std::istringstream mesh_details(grid_info_sample);
        std::istringstream boundary_file(boundary_sample);
        Options opt;

        grid = std::make_shared<CartesianGridTestInterface>(opt,
            std::move(initial_conditions), std::move(mesh_details),
            std::move(boundary_file));
        pf = std::make_shared<PointFunctions>(opt.mach(), opt.gam());
        der = std::make_shared<RegularDerivatives>(
            *pf, grid->shiftX(), grid->shiftY(), grid->dx, grid->dy);
        boundary = std::make_shared<Boundary>(
            *pf, der, opt.reynolds(), opt.prandtl());
        nPointsTotal = grid->nPointsTotal;
        nPointsI = grid->nPointsI;
        nPointsJ = grid->nPointsJ;
    }
    std::shared_ptr<CartesianGridTestInterface> grid;
    std::shared_ptr<RegularDerivatives> der;
    std::shared_ptr<PointFunctions> pf;
    std::shared_ptr<Boundary> boundary;
    int nPointsI, nPointsJ, nPointsTotal;
    double dx, dy;
    Flux lower_x, lower_y;
    Flux upper_x, upper_y;
};

TEST_F(BoundaryTest, TestingMakesSense)
{
    ASSERT_FLOAT_EQ(grid->xmin, 0);
    ASSERT_FLOAT_EQ(grid->dx, 0.25);
    ASSERT_FLOAT_EQ((grid->boundary()).size(), 16);
}

TEST_F(BoundaryTest, ConvectionIsSymmetric)
{
    for (auto& bp : grid->boundary()) {
        int ind = bp.ind;
        int i = grid->indI(ind);
        int j = grid->indJ(ind);
        int new_i = (grid->nPointsI - 1) - i;
        int new_ind = grid->IND(new_i, j);
        for (auto& bp_symmetric : grid->boundary()) {
            auto p = grid->values(bp.ind);
            auto ma = sqrt((pf->ru2(p) + pf->rv2(p)) / pf->rho(p))
                / pf->sound_speed(p);
            ASSERT_LE(ma, 1);
            if (bp_symmetric.ind == new_ind and bp_symmetric.ind > bp.ind) {

                lower_x = boundary->convection_x(*grid, bp, 0.0);
                upper_x = boundary->convection_x(*grid, bp_symmetric, 0.0);
                EXPECT_FLOAT_EQ(lower_x.rho, upper_x.rho);
                EXPECT_FLOAT_EQ(lower_x.ru, upper_x.ru);
                EXPECT_FLOAT_EQ(lower_x.rv, -upper_x.rv);
                EXPECT_FLOAT_EQ(lower_x.e, upper_x.e);
                lower_y = boundary->convection_y(*grid, bp, 0.0);
                upper_y = boundary->convection_y(*grid, bp_symmetric, 0.0);
                EXPECT_FLOAT_EQ(lower_y.rho, upper_y.rho);
                EXPECT_FLOAT_EQ(lower_y.ru, upper_y.ru);
                EXPECT_FLOAT_EQ(lower_y.rv, -upper_y.rv);
                EXPECT_FLOAT_EQ(lower_y.e, upper_y.e);
            }
        }
    }
}

TEST_F(BoundaryTest, DissipationIsSymmetric)
{
    for (auto& bp : grid->boundary()) {
        int ind = bp.ind;
        int i = grid->indI(ind);
        int j = grid->indJ(ind);
        int new_i = (grid->nPointsI - 1) - i;
        int new_ind = grid->IND(new_i, j);
        for (auto& bp_symmetric : grid->boundary()) {
            if (bp_symmetric.ind == new_ind and bp_symmetric.ind > bp.ind) {
                lower_x = boundary->dissipation_x(*grid, bp);
                upper_x = boundary->dissipation_x(*grid, bp_symmetric);
                EXPECT_FLOAT_EQ(lower_x.rho, upper_x.rho);
                EXPECT_FLOAT_EQ(lower_x.ru, upper_x.ru);
                EXPECT_FLOAT_EQ(lower_x.rv, -upper_x.rv);
                EXPECT_FLOAT_EQ(lower_x.e, upper_x.e);

                lower_y = boundary->dissipation_y(*grid, bp);
                upper_y = boundary->dissipation_y(*grid, bp_symmetric);
                EXPECT_FLOAT_EQ(lower_y.rho, upper_y.rho);
                EXPECT_FLOAT_EQ(lower_y.ru, upper_y.ru);
                EXPECT_FLOAT_EQ(lower_y.rv, -upper_y.rv);
                EXPECT_FLOAT_EQ(lower_y.e, upper_y.e);
            }
        }
    }
}
