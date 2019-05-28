#include "../../input_output/options.hpp"
#include "../../utils/operators_overloads.hpp"
#include "../../utils/point_def.hpp"
#include "../../utils/shock_discontinuity_def.hpp"
#include "../../utils/shock_discontinuity_handler.hpp"
#include "../useful_alias.hpp"
#include "gtest/gtest.h"

Options opt;

// Values are picked to be the same as in the shock appearing in sod shock tube
// problem
//

TEST(ShockDiscontinuityTest, PrepareValuesXTest)
{
    ShockHandler sh(opt);
    Point p_left(
        0.265573711705, 0.265573711705 * 0.9274526200489, 0.0, 0.87204449747);

    Point p_right(0.125, 0.0, 0.0, 0.25);
    ShockDiscontinuity shock = sh.create_shock(p_left, p_right, 'x');
    EXPECT_TRUE(shock.is_strong());

    EXPECT_NEAR(-1.75215573203017, shock.w, 1e-5);
}

TEST(ShockDiscontinuityTest, PrepareValuesYTest)
{
    ShockHandler sh(opt);
    Point p_left(
        0.265573711705, 0.0, 0.265573711705 * 0.9274526200489, 0.87204449747);
    Point p_right(0.125, 0.0, 0.0, 0.25);
    ShockDiscontinuity shock = sh.create_shock(p_left, p_right, 'y');
    EXPECT_FALSE(shock.is_weak());
    EXPECT_TRUE(shock.is_strong());
    EXPECT_NEAR(-1.75215573203017, shock.w, 1e-5);
}

TEST(ShockDiscontinuityTest, PrepareValuesXTestReverse)
{
    ShockHandler sh(opt);
    Point p_right(
        0.265573711705, -0.265573711705 * 0.9274526200489, 0.0, 0.87204449747);
    Point p_left(0.125, 0.0, 0.0, 0.25);

    ShockDiscontinuity shock = sh.create_shock(p_left, p_right, 'x');
    EXPECT_TRUE(shock.is_strong());
    EXPECT_NEAR(-1.75215573203017, shock.w, 1e-5);
}

/*Handcrafted test case*/

TEST(ShockDiscontinuityTest, PrepareValuesXHandcraft)
{
    ShockHandler sh(opt);
    Point p_right(1.0, 0.0, 0.0, 25 / 14.);
    Point p_left(8 / 3., 8 / 3. * 5 / 4., 0.0, 425 / 42.);

    ShockDiscontinuity shock = sh.create_shock(p_left, p_right, 'x');
    EXPECT_TRUE(shock.is_strong());
    EXPECT_NEAR(-2, shock.w, 1e-5);
}

TEST(ShockDiscontinuityTest, ExpansionShock)
{
    /*Expand to left*/
    ShockHandler sh(opt);
    Point p_left(
        0.265573711705, -0.265573711705 * 0.9274526200489, 0.0, 0.87204449747);
    Point p_right(0.125, 0.0, 0.0, 0.25);
    ShockDiscontinuity shock = sh.create_shock(p_left, p_right, 'x');
    EXPECT_TRUE(shock.is_weak());

    /*Expand to right*/
    p_right = Point(
        0.265573711705, 0.265573711705 * 0.9274526200489, 0.0, 0.87204449747);
    p_left = Point(0.125, 0.0, 0.0, 0.25);
    shock = sh.create_shock(p_left, p_right, 'x');
    EXPECT_TRUE(shock.is_weak());
}

TEST(ShockDiscontinuityTest, WallCollision)
{
    // Creates a shock
    ShockHandler sh(opt);
    Point p_left(8 / 3., 8 / 3. * 5 / 4., 0.0, 425 / 42.);
    Point p_right(1.0, 0.0, 0.0, 25 / 14.);
    ShockDiscontinuity shock = sh.create_shock(p_left, p_right, 'x');

    sh.compute_wall_interaction(shock, 0.0);
    // Copy to compare before and after prepare
    auto after_col = shock;
    sh.prepare_values(shock);

    EXPECT_TRUE(after_col.left() == shock.left()); // Untouched by prepare
    // Almost the same
    EXPECT_NEAR(after_col.w, shock.w, 1e-7);
    EXPECT_NEAR(after_col.right().rho(), shock.right().rho(), 1e-7);
    EXPECT_NEAR(after_col.right().ru(), shock.right().ru(), 1e-7);
    EXPECT_NEAR(after_col.right().rv(), shock.right().rv(), 1e-7);
    EXPECT_NEAR(after_col.right().e(), shock.right().e(), 1e-7);
}
