#ifndef FLAG_HANDLER_HPP
#define FLAG_HANDLER_HPP
#include <cstdint>

extern const int FLUID_POINT;
extern const int SOLID_POINT;
extern const int GHOST_POINT;
extern const int BOUNDARY_POINT;

namespace flag_functions {
namespace constants {
    extern const int TYPE_MASK_SIZE;
    extern const int TYPE_MASK_SHIFT;
    extern const int TYPE_MASK;

    extern const int DERX_MASK_SIZE;
    extern const int DERX_MASK_SHIFT;
    extern const int DERX_MASK;

    extern const int DERY_MASK_SIZE;
    extern const int DERY_MASK_SHIFT;
    extern const int DERY_MASK;

    extern const int LEFT_MASK;
    extern const int RIGHT_MASK;
}

int point_type(int flag);
int derx(int flag);
int dery(int flag);
int left(int der_flag);
int right(int der_flag);

} /*flag_functions namespace*/

#endif /* FLAG_HANDLER_HPP */
