#include "gtest/gtest.h"

#include "../../input_output/options.hpp"
#include "../../utils/point_functions.hpp"
#include "../regular_derivatives.hpp"

#include <iostream>
#include <memory>
#include <sstream>

#include "../../grid/test/cartesian_grid_test_interface.hpp"
#include "sample_inputs.inc"

/*
 * Grid read from aux file has:
 * xmin = ymin = 0
 * nPointsI = nPointsJ = 11
 * dx = dy = 0.1
 *
 * rho(x,y) = cos( 2x )
 * ru(x,y)  = sin( 3y )
 * rv(x,y)  = exp(x) - exp(y)
 * e(x,y)   = cos( xy )
 */

double drho_dx(double x, double) { return -2 * sin(2 * x); }
double dru_dx(double, double) { return 0.0; }
double drv_dx(double x, double) { return exp(x); }
double de_dx(double x, double y) { return -y * sin(x * y); }

double drho_dy(double, double) { return 0.0; }
double dru_dy(double, double y) { return 3 * cos(3 * y); }
double drv_dy(double, double y) { return -exp(y); }
double de_dy(double x, double y) { return -x * sin(x * y); }

double drho_dxx(double x, double) { return -4 * cos(2 * x); }
double dru_dxx(double, double) { return 0.0; }
double drv_dxx(double x, double) { return exp(x); }
double de_dxx(double x, double y) { return -y * y * cos(x * y); }

double drho_dyy(double, double) { return 0.0; }
double dru_dyy(double, double y) { return -9 * sin(3 * y); }
double drv_dyy(double, double y) { return -exp(y); }
double de_dyy(double x, double y) { return -x * x * cos(x * y); }

double drho_dxy(double, double) { return 0.0; }
double dru_dxy(double, double) { return 0.0; }
double drv_dxy(double, double) { return 0.0; }
double de_dxy(double x, double y) { return -sin(x * y) - x * y * cos(x * y); }

class RegularDerivativesTest : public testing::Test {
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
        nPointsTotal = grid->nPointsTotal;
        nPointsI = grid->nPointsI;
        nPointsJ = grid->nPointsJ;
        dx = grid->dx;
        dy = grid->dy;
    }
    std::shared_ptr<CartesianGridTestInterface> grid;
    std::shared_ptr<RegularDerivatives> der;
    std::shared_ptr<PointFunctions> pf;
    int nPointsI, nPointsJ, nPointsTotal;
    double dx, dy;
};

TEST_F(RegularDerivativesTest, TestingMakesSense)
{
    ASSERT_FLOAT_EQ(grid->xmin, 0);
    ASSERT_FLOAT_EQ(grid->dx, 0.1);
    ASSERT_FLOAT_EQ((grid->boundary()).size(), 40);
}

TEST_F(RegularDerivativesTest, DX)
{
    for (int ind = 0; ind < nPointsTotal; ind += 5) {
        double x = grid->X(ind);
        double y = grid->Y(ind);
        ASSERT_NEAR(der->DX(*grid, alias::RHO, ind), drho_dx(x, y), 2 * dx);
        ASSERT_NEAR(der->DX(*grid, alias::RU, ind), dru_dx(x, y), 2 * dx);
        ASSERT_NEAR(der->DX(*grid, alias::RV, ind), drv_dx(x, y), 2 * dx);
        ASSERT_NEAR(der->DX(*grid, alias::E, ind), de_dx(x, y), 2 * dx);
    }
}

TEST_F(RegularDerivativesTest, DY)
{
    for (int ind = 0; ind < nPointsTotal; ind += 5) {
        double x = grid->X(ind);
        double y = grid->Y(ind);
        ASSERT_NEAR(der->DY(*grid, alias::RHO, ind), drho_dy(x, y), 2 * dy);
        ASSERT_NEAR(der->DY(*grid, alias::RU, ind), dru_dy(x, y), 2 * dy);
        ASSERT_NEAR(der->DY(*grid, alias::RV, ind), drv_dy(x, y), 2 * dy);
        ASSERT_NEAR(der->DY(*grid, alias::E, ind), de_dy(x, y), 2 * dy);
    }
}

TEST_F(RegularDerivativesTest, DXX)
{
    // Bulk
    for (int i = 1; i < nPointsI - 1; i += 3) {
        for (int j = 1; j < nPointsJ - 1; j += 4) {
            int ind = grid->IND(i, j);
            double x = grid->X(ind);
            double y = grid->Y(ind);
            ASSERT_NEAR(
                der->DXX(*grid, alias::RHO, ind), drho_dxx(x, y), 8 * dx * dx);
            ASSERT_NEAR(
                der->DXX(*grid, alias::RU, ind), dru_dxx(x, y), 8 * dx * dx);
            ASSERT_NEAR(
                der->DXX(*grid, alias::RV, ind), drv_dxx(x, y), 8 * dx * dx);
            ASSERT_NEAR(
                der->DXX(*grid, alias::E, ind), de_dxx(x, y), 8 * dx * dx);
        }
    }
    // Borders
    for (int b = 0; b < nPointsI; b++) // Take advantage that grid is square
    {
        int ind_vec[4] = {grid->IND(0, b), grid->IND(nPointsI - 1, b),
            grid->IND(b, 0), grid->IND(b, nPointsJ - 1)};
        for (int i = 0; i < 4; i++) {
            int ind = ind_vec[i];
            double x = grid->X(ind);
            double y = grid->Y(ind);
            ASSERT_NEAR(
                der->DXX(*grid, alias::RHO, ind), drho_dxx(x, y), 8 * dx);
            ASSERT_NEAR(der->DXX(*grid, alias::RU, ind), dru_dxx(x, y), 8 * dx);
            ASSERT_NEAR(der->DXX(*grid, alias::RV, ind), drv_dxx(x, y), 8 * dx);
            ASSERT_NEAR(der->DXX(*grid, alias::E, ind), de_dxx(x, y), 8 * dx);
        }
    }
}

TEST_F(RegularDerivativesTest, DYY)
{
    // Bulk
    for (int i = 1; i < nPointsI - 1; i += 3)
        for (int j = 0; j < nPointsJ - 1; j += 4) {
            int ind = grid->IND(i, j);
            double x = grid->X(ind);
            double y = grid->Y(ind);
            ASSERT_NEAR(
                der->DYY(*grid, alias::RHO, ind), drho_dyy(x, y), 27 * dy * dy);
            ASSERT_NEAR(
                der->DYY(*grid, alias::RU, ind), dru_dyy(x, y), 27 * dy * dy);
            ASSERT_NEAR(
                der->DYY(*grid, alias::RV, ind), drv_dyy(x, y), 27 * dy * dy);
            ASSERT_NEAR(
                der->DYY(*grid, alias::E, ind), de_dyy(x, y), 27 * dy * dy);
        }
    // Borders
    for (int b = 0; b < nPointsI; b++) // Take advantage that grid is square
    {
        int ind_vec[4] = {grid->IND(0, b), grid->IND(nPointsI - 1, b),
            grid->IND(b, 0), grid->IND(b, nPointsJ - 1)};
        for (int i = 0; i < 4; i++) {
            int ind = ind_vec[i];
            double x = grid->X(ind);
            double y = grid->Y(ind);
            ASSERT_NEAR(
                der->DYY(*grid, alias::RHO, ind), drho_dyy(x, y), 27 * dy);
            ASSERT_NEAR(
                der->DYY(*grid, alias::RU, ind), dru_dyy(x, y), 27 * dy);
            ASSERT_NEAR(
                der->DYY(*grid, alias::RV, ind), drv_dyy(x, y), 27 * dy);
            ASSERT_NEAR(der->DYY(*grid, alias::E, ind), de_dyy(x, y), 27 * dy);
        }
    }
}

TEST_F(RegularDerivativesTest, DXY)
{
    for (int ind = 0; ind < nPointsTotal; ind += 5) {
        double x = grid->X(ind);
        double y = grid->Y(ind);
        ASSERT_NEAR(der->DXY(*grid, alias::RHO, ind), drho_dxy(x, y), 2 * dy);
        ASSERT_NEAR(der->DXY(*grid, alias::RU, ind), dru_dxy(x, y), 2 * dy);
        ASSERT_NEAR(der->DXY(*grid, alias::RV, ind), drv_dxy(x, y), 2 * dy);
        ASSERT_NEAR(der->DXY(*grid, alias::E, ind), de_dxy(x, y), 2 * dy);
    }
}
