#include "gtest/gtest.h"

#include <iostream>
#include <memory>
#include <sstream>

#include "../../boundary/boundary.hpp"
#include "../../derivatives/regular_derivatives.hpp"
#include "../../grid/test/cartesian_grid_test_interface.hpp"
#include "sample_inputs.inc"

class TimeIntegratorTest : public testing::Test {
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
        boundary = std::make_shared<Boundary>(*pf, *der);
        nPointsTotal = grid->nPointsTotal;
        nPointsI = grid->nPointsI;
        nPointsJ = grid->nPointsJ;
    }
    std::shared_ptr<CartesianGridTestInterface> grid;
    std::shared_ptr<RegularDerivatives> der;
    std::shared_ptr<PointFunctions> pf;
    std::shared_ptr<Boundary> boundary;
    int nPointsI, nPointsJ, nPointsTotal;
    Flux result;
};
