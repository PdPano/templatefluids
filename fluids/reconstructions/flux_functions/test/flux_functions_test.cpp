#include "../"
#include "gtest/gtest.h"

TEST(FlagHandlerTest, Basic)
{
    uint16_t v1 = 1, v2 = 2, v3 = 3, v4 = 4;
    using namespace flag_functions;
    uint16_t sample_flag
        = (v1 << constants::TYPE_MASK_SHIFT) // type
        + ((v1 << constants::DERX_MASK_SHIFT)
              + (v2 << (constants::DERX_MASK_SHIFT + 3))) // DerX
        + ((v3 << constants::DERY_MASK_SHIFT)
              + (v4 << (constants::DERY_MASK_SHIFT + 3))); // DerY

    EXPECT_EQ(point_type(sample_flag), v1);
    EXPECT_EQ(left(derx(sample_flag)), v1);
    EXPECT_EQ(right(derx(sample_flag)), v2);
    EXPECT_EQ(left(dery(sample_flag)), v3);
    EXPECT_EQ(right(dery(sample_flag)), v4);
}
