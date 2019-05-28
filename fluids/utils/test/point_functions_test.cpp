#include "../point_functions.hpp"
#include "../../utils/point_def.hpp"
#include "../useful_alias.hpp"
#include "gtest/gtest.h"
#include <cmath>

PointFunctions p_funcs(5.0, 1.4);

Point p(1.0, -1.0, 2.0, 3.5);

// Values are cherry-picked to be pretty:
//
// rho_internal_e = E - (ru^2 + rv^2)/2
//           = 3.5 - 2.5 = 1.0
//
// p = (gam-1) rho_internal_e = 0.4
//
// T = gam*mach^2*p/rho = 1.4*5^2*0.4/1 = 14
//
// c = sqrt( gam*p/rho ) = sqrt(1.4*0.4) = sqrt(0.56)

TEST(PointFunctionsTest, TrivialFunctions)
{
    EXPECT_FLOAT_EQ(p_funcs.rho(p), 1.0);
    EXPECT_FLOAT_EQ(p_funcs.ru(p), -1.0);
    EXPECT_FLOAT_EQ(p_funcs.rv(p), 2.0);
    EXPECT_FLOAT_EQ(p_funcs.e(p), 3.5);
}

TEST(PointFunctionsTest, BasicFunctions)
{
    EXPECT_FLOAT_EQ(p_funcs.u(p), -1.0);
    EXPECT_FLOAT_EQ(p_funcs.v(p), 2.0);
    EXPECT_FLOAT_EQ(p_funcs.ru2(p), 1.0);
    EXPECT_FLOAT_EQ(p_funcs.rv2(p), 4.0);
}

TEST(PointFunctionsTest, ThermodynamicalFunctions)
{
    EXPECT_FLOAT_EQ(p_funcs.temperature(p), 14.0);
    EXPECT_FLOAT_EQ(p_funcs.pressure(p), 0.4);
    EXPECT_FLOAT_EQ(p_funcs.sound_speed(p), sqrt(0.56));
    EXPECT_FLOAT_EQ(p_funcs.e_from_T(p, 14), 3.5);
}

TEST(PointFunctionsAliasTest, Access)
{
    EXPECT_FLOAT_EQ((p_funcs.*alias::RHO)(p), 1.0);
}
