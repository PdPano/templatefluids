#include "flag_handler.hpp"

const int FLUID_POINT = 0;
const int SOLID_POINT = 1;
const int GHOST_POINT = 2;
const int BOUNDARY_POINT = 3;

namespace flag_functions {
namespace constants {
    const int TYPE_MASK_SIZE = 2;
    const int TYPE_MASK_SHIFT = 0;
    const int TYPE_MASK = ((1 << TYPE_MASK_SIZE) - 1) << TYPE_MASK_SHIFT;

    const int DERX_MASK_SIZE = 8;
    const int DERX_MASK_SHIFT = TYPE_MASK_SIZE + TYPE_MASK_SHIFT;
    const int DERX_MASK = ((1 << DERX_MASK_SIZE) - 1) << DERX_MASK_SHIFT;

    const int DERY_MASK_SIZE = 8;
    const int DERY_MASK_SHIFT = DERX_MASK_SIZE + DERX_MASK_SHIFT;
    const int DERY_MASK = ((1 << DERY_MASK_SIZE) - 1) << DERY_MASK_SHIFT;

    const int LEFT_MASK = 15;
    const int RIGHT_MASK = (15 << 4);
} // namespace constants

int point_type(int flag)
{
    return (flag & constants::TYPE_MASK) >> constants::TYPE_MASK_SHIFT;
}

int derx(int flag)
{
    return (flag & constants::DERX_MASK) >> constants::DERX_MASK_SHIFT;
}

int dery(int flag)
{
    return (flag & constants::DERY_MASK) >> constants::DERY_MASK_SHIFT;
}

int left(int der_flag) { return der_flag & constants::LEFT_MASK; }

int right(int der_flag) { return (der_flag & constants::RIGHT_MASK) >> 4; }
} // namespace flag_functions
