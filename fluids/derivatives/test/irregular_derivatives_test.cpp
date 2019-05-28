#include "gtest/gtest.h"

#include "../../grid/karagiozis_grid.hpp"
#include "../../input_output/options.hpp"
#include "../../input_output/readers/reader.hpp"
#include "../irregular_derivatives.hpp"

#include <iostream>
#include <memory>
#include <sstream>

#include "sample_inputs_irregular.inc"

/*
 * Grid read from aux file has:
 * xmin = ymin = 0
 * nPointsI = nPointsJ = 5
 * dx = dy = 0.25
 *
 * rho(x,y) = cos( x ) * sin( y ) + 2
 * u = v = 0
 * e = 20
 */

double drho_dx(double x, double y) { return -sin(x) * sin(y); }

double drho_dy(double x, double y) { return cos(x) * cos(y); }

double drho_dxx(double x, double y) { return -cos(x) * sin(y); }

double drho_dyy(double x, double y) { return -cos(x) * sin(y); }

double drho_dxy(double x, double y) { return -sin(x) * cos(y); }

class IrregularDerivativesTest : public testing::Test {
protected:
    virtual void SetUp()
    {

        std::istringstream initial_conditions(initial_conditions_sample);
        std::istringstream mesh_details(grid_info_sample);
        std::istringstream boundary_file(boundary_sample);
        std::istringstream immersed_file(immersed_interface_sample);
        Options opt;

        Reader reader(opt, std::move(initial_conditions),
            std::move(mesh_details), std::move(boundary_file),
            std::move(immersed_file));
        grid = std::make_shared<KaragiozisGrid>(reader, opt);
        grid->grid_specific_update();
        pf = std::make_shared<PointFunctions>(opt.mach(), opt.gam());
        der = std::make_shared<IrregularDerivatives>(
            *pf, grid->shiftX(), grid->shiftY(), grid->dx, grid->dy);
        nPointsTotal = grid->nPointsTotal;
        nPointsI = grid->nPointsI;
        nPointsJ = grid->nPointsJ;
        dx = grid->dx;
        dy = grid->dy;
    }
    std::shared_ptr<KaragiozisGrid> grid;
    std::shared_ptr<IrregularDerivatives> der;
    std::shared_ptr<PointFunctions> pf;
    int nPointsI, nPointsJ, nPointsTotal;
    double dx, dy;
    double x, y;
};

TEST_F(IrregularDerivativesTest, TestingMakesSense)
{
    ASSERT_FLOAT_EQ(grid->xmin, 0);
    ASSERT_FLOAT_EQ(grid->dx, 0.25);
}

TEST_F(IrregularDerivativesTest, DX)
{
    auto D = [&](int ind) { return der->DX(*grid, alias::RHO, ind); };
    auto Exact = [&](int ind) {
        return drho_dx(dynamic_cast<CartesianGrid*>(grid.get())->X(ind),
            dynamic_cast<CartesianGrid*>(grid.get())->Y(ind));
    };
    EXPECT_EQ(D(5), 0.0);            // Left Boundary with disc on the right
    EXPECT_NEAR(D(6), Exact(6), dx); // One disc to the left (far)
    EXPECT_NEAR(D(8), Exact(8), dx); // One disc to the right (close)
    EXPECT_EQ(D(9), 0.0);            // Right Boundary with disc on the left

    auto disc = grid->first_disc_x(12);
    disc->left().set_rho(1.0);
    disc->right().set_rho(1.0);
    EXPECT_EQ(D(12), 0.0); // Discontinuities on either side, too close
    EXPECT_NE(D(13), 0.0); // Discontinuities on either side, not too close
}

TEST_F(IrregularDerivativesTest, DY)
{
    auto D = [&](int ind) { return der->DY(*grid, alias::RHO, ind); };
    auto Exact = [&](int ind) {
        return drho_dy(dynamic_cast<CartesianGrid*>(grid.get())->X(ind),
            dynamic_cast<CartesianGrid*>(grid.get())->Y(ind));
    };
    EXPECT_EQ(D(1), 0.0);                // Left Boundary with disc on the right
    EXPECT_NEAR(D(6), Exact(6), 2 * dy); // One disc to the left (far)
    EXPECT_NEAR(D(16), Exact(16), dy);   // One disc to the right (close)
    EXPECT_EQ(D(21), 0.0);               // Right Boundary with disc on the left

    auto disc = grid->first_disc_y(12);
    disc->left().set_rho(1.0);
    disc->right().set_rho(1.0);
    EXPECT_EQ(D(12), 0.0); // Discontinuities on either side, too close
    EXPECT_NE(D(17), 0.0); // Discontinuities on either side, not too close
}

TEST_F(IrregularDerivativesTest, DXX)
{
    auto D = [&](int ind) { return der->DXX(*grid, alias::RHO, ind); };
    auto Exact = [&](int ind) {
        return drho_dxx(dynamic_cast<CartesianGrid*>(grid.get())->X(ind),
            dynamic_cast<CartesianGrid*>(grid.get())->Y(ind));
    };

    EXPECT_EQ(D(5), 0.0);            // Left Boundary with disc on the right
    EXPECT_NEAR(D(6), Exact(6), dx); // One disc to the left (far)
    EXPECT_NEAR(D(7), Exact(7), dx * dx); // Should use regular
    EXPECT_NEAR(D(8), Exact(8), dx);      // One disc to the right (close)
    EXPECT_EQ(D(9), 0.0); // Right boundary with disc on the left

    // Left Boundary with disc on the second
    auto disc0 = grid->first_disc_x(11);
    disc0->left().set_rho(2.42634964);       // Real value
    EXPECT_LE(fabs(D(10)), fabs(Exact(10))); // Underestimates
    EXPECT_EQ(D(10), D(11));                 // Two points with same der

    EXPECT_EQ(D(12), 0.0); // Discs on both sides

    // Between two discs
    auto disc1 = grid->first_disc_x(16);
    disc1->right().set_rho(2.606176388); // Real value
    disc1 = grid->first_disc_x(18);
    disc1->left().set_rho(2.4869771);        // Real value
    EXPECT_EQ(D(17), D(18));                 // Two equal
    EXPECT_LE(fabs(D(18)), fabs(Exact(18))); // Underestimate
    EXPECT_GE(D(18) * Exact(18), 0);         // Same sign

    EXPECT_EQ(D(22), D(23));          // Between two discs
    EXPECT_NEAR(D(22), 0.0, dx * dx); // Linear interpolation
}

TEST_F(IrregularDerivativesTest, DYY)
{
    auto D = [&](int ind) { return der->DYY(*grid, alias::RHO, ind); };
    auto Exact = [&](int ind) {
        return drho_dyy(dynamic_cast<CartesianGrid*>(grid.get())->X(ind),
            dynamic_cast<CartesianGrid*>(grid.get())->Y(ind));
    };

    EXPECT_EQ(D(1), 0.0);                // Left Boundary with disc on the right
    EXPECT_NEAR(D(6), Exact(6), 5 * dy); // One disc to the left (far)
    EXPECT_NEAR(D(11), Exact(11), dy * dy); // Should use regular
    EXPECT_NEAR(D(16), Exact(16), dy);      // One disc to the right (close)
    EXPECT_EQ(D(21), 0.0); // Right boundary with disc on the left

    // Left Boundary with disc on the second
    auto disc0 = grid->first_disc_y(7);
    disc0->left().set_rho(2.40135224612);
    EXPECT_LE(fabs(D(7)), fabs(Exact(7))); // Underestimates
    EXPECT_EQ(D(2), D(7));                 // Two points with same der

    EXPECT_EQ(D(12), 0.0); // Discs on both sides

    // Between two discs
    auto disc1 = grid->first_disc_y(8);
    disc1->right().set_rho(2.334629); // Real value
    disc1 = grid->first_disc_y(18);
    disc1->left().set_rho(2.5119744);        // Real value
    EXPECT_EQ(D(13), D(18));                 // Two equal
    EXPECT_LE(fabs(D(18)), fabs(Exact(18))); // Underestimate
    EXPECT_GE(D(18) * Exact(18), 0);         // Same sign

    EXPECT_EQ(D(14), D(19));          // Between two discs
    EXPECT_NEAR(D(14), 0.0, dy * dy); // Linear interpolation
}

TEST(IrregularDerivativesTestCross, DXY)
{
    std::istringstream initial_conditions(initial_conditions_sample);
    std::istringstream mesh_details(grid_info_sample);
    std::istringstream boundary_file(boundary_sample);
    std::istringstream immersed_file(immersed_interface_cross);
    Options opt;

    Reader reader(opt, std::move(initial_conditions), std::move(mesh_details),
        std::move(boundary_file), std::move(immersed_file));
    KaragiozisGrid grid(reader, opt);
    double dx = grid.dx;
    double dy = grid.dy;
    PointFunctions pf(opt.mach(), opt.gam());
    IrregularDerivatives der(pf, grid.shiftX(), grid.shiftY(), dx, dy);
    grid.grid_specific_update();

    auto D = [&](int ind) { return der.DXY(grid, alias::RHO, ind); };
    auto Exact = [&](int ind) {
        return drho_dxy(dynamic_cast<CartesianGrid*>(&grid)->X(ind),
            dynamic_cast<CartesianGrid*>(&grid)->Y(ind));
    };

    grid.setRho(1.0, 0);
    EXPECT_NEAR(D(6), Exact(6), dx + dy);

    grid.setRho(1.0, 20);
    EXPECT_NEAR(D(16), Exact(16), dx + dy);

    EXPECT_EQ(D(13), 0.0);
    EXPECT_NEAR(D(18), Exact(18), dx + dy);
}
