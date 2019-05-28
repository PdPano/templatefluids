#include "../ghias_grid.hpp"
#include "../../input_output/options.hpp"
#include "../../input_output/readers/reader.hpp"
#include "../../input_output/stream_from_file.hpp"
#include "gtest/gtest.h"

#include <sstream>

#include "sample_inputs.inc"

TEST(GhiasGridTest, testInterpolationAllFluid)
{
    std::istringstream initial_conditions(initial_conditions_sample);
    std::istringstream mesh_details(grid_info_sample);
    std::istringstream boundary_file(boundary_sample);
    std::istringstream immersed_file(immersed_interface_all_fluid);
    Options opt;

    Reader reader(opt, std::move(initial_conditions), std::move(mesh_details),
        std::move(boundary_file), std::move(immersed_file));
    GhiasGrid grid(reader, opt);

    grid.grid_specific_update();

    EXPECT_FLOAT_EQ(grid.rho(2), 1.5);
    EXPECT_FLOAT_EQ(grid.ru(2), 0);
    EXPECT_FLOAT_EQ(grid.rv(2), 0);
    EXPECT_FLOAT_EQ(grid.e(2), 45);
}

TEST(GhiasGridTest, testInterpolationOneBorder)
{
    std::istringstream initial_conditions(initial_conditions_sample);
    std::istringstream mesh_details(grid_info_sample);
    std::istringstream boundary_file(boundary_sample);
    std::istringstream immersed_file(immersed_interface_one_border);
    Options opt;

    Reader reader(opt, std::move(initial_conditions), std::move(mesh_details),
        std::move(boundary_file), std::move(immersed_file));
    GhiasGrid grid(reader, opt);

    grid.grid_specific_update();

    EXPECT_FLOAT_EQ(grid.rho(2), 1.333333);
    EXPECT_FLOAT_EQ(grid.ru(2), 0);
    EXPECT_FLOAT_EQ(grid.rv(2), 0);
    EXPECT_FLOAT_EQ(grid.e(2), 44.444444);
}
