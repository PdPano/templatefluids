#include "../karagiozis_grid.hpp"
#include "../../input_output/options.hpp"
#include "../../input_output/readers/reader.hpp"
#include "../../input_output/stream_from_file.hpp"
#include "gtest/gtest.h"

#include <sstream>
#include <tuple>

#include "sample_inputs_karagiozis.inc"

double analytical_rho(double x, double y) { return cos(x) * sin(y) + 2; }

auto rho_left_and_right(GenericDiscontinuity* bd)
{
    return std::make_tuple(bd->left().rho(), bd->right().rho());
}

TEST(KaragiozisGridTest, testMultipleDiscontinuities)
{
    std::istringstream initial_conditions(initial_conditions_sample);
    std::istringstream mesh_details(grid_info_sample);
    std::istringstream boundary_file(boundary_sample);
    std::istringstream immersed_file(
        immersed_interface_multiple_discontinuities);
    Options opt;

    Reader reader(opt, std::move(initial_conditions), std::move(mesh_details),
        std::move(boundary_file), std::move(immersed_file));
    KaragiozisGrid grid(reader, opt);

    auto xmap = grid.discontinuity_map_x();
    EXPECT_EQ(xmap->find(0), xmap->end());
    EXPECT_EQ((*xmap)[7].size(), 3);
    EXPECT_EQ((*xmap)[7][0]->frac, 0.2);
    EXPECT_EQ((*xmap)[7][1]->frac, 0.5);
    EXPECT_EQ((*xmap)[7][2]->frac, 0.6);

    auto ymap = grid.discontinuity_map_y();
    EXPECT_EQ(ymap->find(0), ymap->end());
    EXPECT_EQ((*ymap)[7].size(), 3);
    EXPECT_EQ((*ymap)[7][0]->frac, 0.2);
    EXPECT_EQ((*ymap)[7][1]->frac, 0.5);
    EXPECT_EQ((*ymap)[7][2]->frac, 0.7);
}

TEST(KaragiozisGridTest, testExtrapolation)
{
    std::istringstream initial_conditions(initial_conditions_sample);
    std::istringstream mesh_details(grid_info_sample);
    std::istringstream boundary_file(boundary_sample);
    std::istringstream immersed_file(immersed_interface_extrapolation);
    Options opt;

    Reader reader(opt, std::move(initial_conditions), std::move(mesh_details),
        std::move(boundary_file), std::move(immersed_file));
    KaragiozisGrid grid(reader, opt);

    const double dx = grid.dx;
    const double dy = grid.dy;

    grid.grid_specific_update();

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

    EXPECT_EQ(grid.to_revisit().size(), 9);
}
