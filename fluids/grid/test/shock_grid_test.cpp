#include "../shock_grid.hpp"
#include "../../utils/debug_def.hpp"
#include "../../input_output/options.hpp"
#include "../../input_output/readers/reader.hpp"
#include "../../input_output/stream_from_file.hpp"
#include "../../utils/operators_overloads.hpp"
#include "gtest/gtest.h"

#include <sstream>
#include <tuple>

#include "sample_inputs_shock.inc"

double analytical_rho(double x, double y) { return cos(x) * sin(y) + 2; }

auto rho_left_and_right(GenericDiscontinuity* bd)
{
    return std::make_tuple(bd->left().rho(), bd->right().rho());
}

TEST(ShockGridTest, testSimultaneousReading)
{
    std::istringstream initial_conditions(initial_conditions_sample);
    std::istringstream mesh_details(grid_info_sample);
    std::istringstream boundary_file(boundary_empty);
    std::istringstream immersed_file(immersed_interface_simple);
    std::istringstream shock_file(shock_points_simple);
    Options opt;

    Reader reader(opt, std::move(initial_conditions), std::move(mesh_details),
        std::move(boundary_file), std::move(immersed_file),
        std::move(shock_file));
    ShockGrid grid(reader, opt);

    auto xmap = grid.discontinuity_map_x();
    auto ymap = grid.discontinuity_map_y();

    EXPECT_EQ(xmap->size(), 4);
    EXPECT_EQ(ymap->size(), 2);
}

TEST(ShockGridTest, testShockAngles)
{
    std::istringstream initial_conditions(initial_conditions_sample);
    std::istringstream mesh_details(grid_info_sample);
    std::istringstream boundary_file(boundary_empty);
    std::istringstream immersed_file(immersed_interface_simple);
    std::istringstream shock_file(shock_angles);
    Options opt;

    Reader reader(opt, std::move(initial_conditions), std::move(mesh_details),
        std::move(boundary_file), std::move(immersed_file),
        std::move(shock_file));
    ShockGrid grid(reader, opt);

    grid.compute_shock_angles();

    auto& shocks = grid.shock_points();

    // Indices are in the same order as read
    EXPECT_EQ(shocks[0].ind, 10);
    EXPECT_EQ(shocks[1].ind, 11);
    EXPECT_EQ(shocks[2].ind, 12);
    EXPECT_EQ(shocks[3].ind, 17);
    EXPECT_EQ(shocks[4].ind, 18);
    EXPECT_EQ(shocks[5].ind, 23);

    // All but the first are connected
    EXPECT_FALSE(shocks[0].is_connected);
    for (int i = 1; i < 6; i++) {
        EXPECT_TRUE(shocks[i].is_connected);
    }

    // Three are connected in the same line
    EXPECT_NEAR(shocks[1].theta, -M_PI / 4, 1e-5);
    EXPECT_NEAR(shocks[2].theta, -M_PI / 4, 1e-5);
    EXPECT_NEAR(shocks[3].theta, -M_PI / 4, 1e-5);

    // Two are connected in a different line
    auto real_angle = atan(2.) + M_PI / 2;
    EXPECT_NEAR(shocks[4].theta, real_angle, 1e-5);
    EXPECT_NEAR(shocks[5].theta, real_angle, 1e-5);
}

TEST(ShockGridTest, testShockExtrapolation)
{
    std::istringstream initial_conditions(initial_conditions_sample);
    std::istringstream mesh_details(grid_info_sample);
    std::istringstream boundary_file(boundary_empty);
    std::istringstream immersed_file(std::string("KARAGIOZIS\n0"));
    std::istringstream shock_file(shock_extrapolation);
    Options opt;

    Reader reader(opt, std::move(initial_conditions), std::move(mesh_details),
        std::move(boundary_file), std::move(immersed_file),
        std::move(shock_file));
    ShockGrid grid(reader, opt);

    const double dx = grid.dx;
    const double dy = grid.dy;

    grid.extrapolate_to_shocks();

    double expected_rho;
    double rho_discont_left;
    double rho_discont_right;
    auto xmap = grid.discontinuity_map_x();
    auto ymap = grid.discontinuity_map_y();

    GenericDiscontinuity* discont;

    discont = (*xmap)[17][0];
    expected_rho = analytical_rho(grid.disc_X(discont), grid.disc_Y(discont));
    std::tie(rho_discont_left, rho_discont_right) = rho_left_and_right(discont);
    EXPECT_NEAR(rho_discont_left, expected_rho, dx * dx);
    EXPECT_EQ(rho_discont_right, 3.0);

    discont = (*xmap)[17].back();
    std::tie(rho_discont_left, rho_discont_right) = rho_left_and_right(discont);
    EXPECT_EQ(rho_discont_left, 1.0);
    EXPECT_EQ(rho_discont_right, grid.rho(18));

    discont = (*xmap)[18][0];
    std::tie(rho_discont_left, rho_discont_right) = rho_left_and_right(discont);
    EXPECT_EQ(rho_discont_left, grid.rho(18));
    EXPECT_EQ(rho_discont_right, grid.rho(19));

    discont = (*ymap)[18].back();
    expected_rho = analytical_rho(grid.disc_X(discont), grid.disc_Y(discont));
    std::tie(rho_discont_left, rho_discont_right) = rho_left_and_right(discont);
    EXPECT_NEAR(rho_discont_left, expected_rho, dy * dy);
    EXPECT_EQ(rho_discont_right, grid.rho(23));
}

TEST(ShockGridTest, testShockMovement)
{
    std::istringstream initial_conditions(initial_conditions_sample);
    std::istringstream mesh_details(grid_info_sample);
    std::istringstream boundary_file(boundary_outlet_sample);
    std::istringstream immersed_file(empty_interface);
    std::istringstream shock_file(shock_movement);
    Options opt;

    Reader reader(opt, std::move(initial_conditions), std::move(mesh_details),
        std::move(boundary_file), std::move(immersed_file),
        std::move(shock_file));

    ShockGrid grid(reader, opt);

    auto& spts = grid.shock_points();

    spts[0].w = -1.0;
    spts[1].w = -1.0;
    spts[2].w = 1.0;
    spts[3].w = 1.0;
    spts[4].w = 1.0;
    spts[5].w = 1.0;

    grid.move_shocks(grid.dx);
    std::vector<int> check_x_dir, check_y_dir;
    grid.update_crossed_grid_points(check_x_dir, check_y_dir);

    // Where should check for connecting shocks
    ASSERT_EQ(check_x_dir.size(), 3u);
    EXPECT_EQ(check_x_dir[0], 1);
    EXPECT_EQ(check_x_dir[1], 15);
    EXPECT_EQ(check_x_dir[2], 23);

    ASSERT_EQ(check_y_dir.size(), 3u);
    EXPECT_EQ(check_y_dir[0], 5);
    EXPECT_EQ(check_y_dir[1], 12);
    EXPECT_EQ(check_y_dir[2], 19);

    EXPECT_FALSE(spts[0].is_connected); // Exits bottom
    EXPECT_FALSE(spts[1].is_connected); // Exits left
    EXPECT_TRUE(spts[2].is_connected);  // Still in domain
    EXPECT_TRUE(spts[3].is_connected);  // Still in domain
    EXPECT_FALSE(spts[4].is_connected); // Exits right
    EXPECT_FALSE(spts[5].is_connected); // Exits top

    EXPECT_EQ(spts[0].right(), grid.values(5));
    EXPECT_EQ(spts[1].right(), grid.values(1));
    EXPECT_EQ(spts[2].left(), grid.values(15)); // Cross near boundary
    EXPECT_EQ(spts[3].left(), grid.values(12)); // Cross in bulk
    EXPECT_EQ(spts[4].left(), grid.values(19));
    EXPECT_EQ(spts[5].left(), grid.values(23));
}

TEST(ShockGridTest, testShockMerging)
{
    std::istringstream initial_conditions(initial_conditions_sample);
    std::istringstream mesh_details(grid_info_sample);
    std::istringstream boundary_file(boundary_outlet_sample);
    std::istringstream immersed_file(empty_interface);
    std::istringstream shock_file(shock_merge);
    Options opt;

    Reader reader(opt, std::move(initial_conditions), std::move(mesh_details),
        std::move(boundary_file), std::move(immersed_file),
        std::move(shock_file));
    ShockGrid grid(reader, opt);
    auto* map_x = grid.discontinuity_map_x();

    grid.merge_compatible_shocks();
    bool removed_shocks = grid.remove_disconnected_shocks();
    EXPECT_TRUE(removed_shocks);
    grid.clear_discontinuity_map();
    grid.fill_discontinuity_map();

    EXPECT_EQ((*map_x)[0].size(), 1);
    EXPECT_FLOAT_EQ((*map_x)[0][0]->frac, 0.5);
    EXPECT_EQ((*map_x)[1].size(), 2);
    EXPECT_EQ((*map_x)[2].size(), 1);
    EXPECT_FLOAT_EQ((*map_x)[2][0]->frac, 0.5);
    EXPECT_EQ((*map_x)[5].size(), 2);
    EXPECT_FLOAT_EQ((*map_x)[5][0]->frac, 0.4);
    EXPECT_EQ((*map_x)[6].size(), 2);
    EXPECT_FLOAT_EQ((*map_x)[6][0]->frac, 0.2);
    EXPECT_FLOAT_EQ((*map_x)[6][1]->frac, 0.3);
    EXPECT_EQ((*map_x)[7].size(), 2);
    EXPECT_FLOAT_EQ((*map_x)[7][0]->frac, 0.3);
}

TEST(ShockGridTest, testRemoveShockNearWall)
{
    std::istringstream initial_conditions(initial_conditions_sample);
    std::istringstream mesh_details(grid_info_sample);
    std::istringstream boundary_file(boundary_empty);
    std::istringstream immersed_file(immersed_interface_simple);
    std::istringstream shock_file(shock_points_near_wall);
    Options opt;

    Reader reader(opt, std::move(initial_conditions), std::move(mesh_details),
        std::move(boundary_file), std::move(immersed_file),
        std::move(shock_file));
    ShockGrid grid(reader, opt);

    EXPECT_EQ(grid.shock_points().size(), 4);

    grid.delete_shock_near_wall();
    grid.remove_disconnected_shocks();

    EXPECT_EQ(grid.shock_points().size(), 2);
}
