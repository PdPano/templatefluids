#include "../../input_output/options.hpp"
#include "../../input_output/readers/reader.hpp"
#include "../../input_output/stream_from_file.hpp"
#include "cartesian_grid_test_interface.hpp"
#include "gtest/gtest.h"

#include <sstream>

#include "sample_inputs.inc"

TEST(CartesianGridTest, testDefaultConstructor)
{
    CartesianGridTestInterface grid;
    int ind = 12;
    double v[4] = {1.0, 1.1, 1.2, 1.3};

    grid.setRho(v[0], ind);
    grid.setRU(v[1], ind);
    grid.setRV(v[2], ind);
    grid.setE(v[3], ind);

    ASSERT_EQ(grid.rho(ind), v[0]);
    ASSERT_EQ(grid.ru(ind), v[1]);
    ASSERT_EQ(grid.rv(ind), v[2]);
    ASSERT_EQ(grid.e(ind), v[3]);
}

TEST(CartesianGridTest, testReaderConstructor)
{
    std::istringstream initial_conditions(initial_conditions_sample);
    std::istringstream mesh_details(grid_info_sample);
    std::istringstream boundary_file(boundary_sample);
    Options opt;

    CartesianGridTestInterface grid(opt, std::move(initial_conditions),
        std::move(mesh_details), std::move(boundary_file));

    // Reading grid values correctly
    ASSERT_EQ(grid.rho(0), 1.0);
    ASSERT_EQ(grid.ru(0), 0.0);
    ASSERT_EQ(grid.rv(0), 0.0);
    ASSERT_EQ(grid.e(0), 40.0);

    ASSERT_EQ(grid.rho(1), 2.0);
    ASSERT_EQ(grid.ru(1), 0.0);
    ASSERT_EQ(grid.rv(1), 0.0);
    ASSERT_EQ(grid.e(1), 40.0);

    // Reading flags correctly
    ASSERT_EQ(grid.flag(0), 0);
    ASSERT_EQ(grid.flag(1), 1);

    // Initializing constants correctly
    ASSERT_EQ(grid.nPointsI, 2);
    ASSERT_EQ(grid.nPointsJ, 2);

    // Reading boundary correctly
    auto& boundary = grid.boundary();
    ASSERT_EQ(boundary.size(), 1);
    ASSERT_EQ(boundary[0].v, 3.0);
}
