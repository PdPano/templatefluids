#include "../../../input_output/options.hpp"
#include "../luisa_detector_23.hpp"
#include "../luisa_detector_345.hpp"
#include "../../../grid/test/cartesian_grid_test_interface.hpp"
#include "gtest/gtest.h"
#include <string>

#include "input_sample.inc"

#include <iostream>
#include <cstdio>
#include <fstream>

Options opt;

class LuisaDetectorTest : public testing::Test {
protected:
    virtual void SetUp()
    {
        std::istringstream initial_conditions(initial_conditions_sample);
        std::istringstream mesh_details(grid_info_sample);
        std::istringstream boundary_file(boundary_sample);
        grid = std::make_shared<CartesianGridTestInterface>(opt,
            std::move(initial_conditions), std::move(mesh_details),
            std::move(boundary_file));
    }
    std::shared_ptr<CartesianGridTestInterface> grid;
};

TEST_F(LuisaDetectorTest, DetectsShocks23)
{
    LuisaDetector23 detector(0.5, grid->nPointsI, grid->nPointsJ);
    detector.detect_shocks(*(grid.get()));
    auto& points = detector.set_of_shocked_points();

    EXPECT_TRUE(points.find(9) != points.end());
    EXPECT_TRUE(points.find(10) != points.end());

    EXPECT_TRUE(points.find(279) != points.end());
    EXPECT_TRUE(points.find(310) != points.end());

    printf(" ");
    for (int j = 0; j < grid->nPointsJ; j++) {
        printf(" %02d", j);
    }
    printf("\n");
    for (int i = grid->nPointsI - 1; i >= 0; i--) {
        printf("%02d", i);
        for (int j = 0; j < grid->nPointsJ; j++) {
            auto ind = grid->IND(i, j);
            if (points.find(ind) == points.end()) {
                std::cout << " . ";
            }
            else {
                std::cout << " X ";
            }
        }
        std::cout << std::endl;
    }

    for (auto& p : points) {
        std::cout << " " << p;
    }
    std::cout << std::endl;
}

TEST_F(LuisaDetectorTest, DetectsShocks345)
{
    LuisaDetector345 detector(1.5, grid->nPointsI, grid->nPointsJ);
    detector.detect_shocks(*(grid.get()));
    auto& points = detector.set_of_shocked_points();

    EXPECT_TRUE(points.find(9) != points.end());
    EXPECT_TRUE(points.find(10) != points.end());

    EXPECT_TRUE(points.find(279) != points.end());
    EXPECT_TRUE(points.find(310) != points.end());

    printf(" ");
    for (int j = 0; j < grid->nPointsJ; j++) {
        printf(" %02d", j);
    }
    printf("\n");
    for (int i = grid->nPointsI - 1; i >= 0; i--) {
        printf("%02d", i);
        for (int j = 0; j < grid->nPointsJ; j++) {
            auto ind = grid->IND(i, j);
            if (points.find(ind) == points.end()) {
                std::cout << " . ";
            }
            else {
                std::cout << " X ";
            }
        }
        std::cout << std::endl;
    }

    for (auto& p : points) {
        std::cout << " " << p;
    }
    std::cout << std::endl;
    /* for (auto& p : points) {
         DUMP(p);
     }*/
}
