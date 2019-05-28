#include "useful_alias.hpp"
#include "point_functions.hpp"
namespace alias {

PointProperty RHO = &PointFunctions::rho;
PointProperty RU = &PointFunctions::ru;
PointProperty RV = &PointFunctions::rv;
PointProperty E = &PointFunctions::e;
PointProperty P = &PointFunctions::pressure;
PointProperty C = &PointFunctions::sound_speed;
PointProperty T = &PointFunctions::temperature;
PointProperty RUV = &PointFunctions::ruv;
PointProperty U = &PointFunctions::u;
PointProperty V = &PointFunctions::v;
PointProperty RU2 = &PointFunctions::ru2;
PointProperty RV2 = &PointFunctions::rv2;
PointProperty INTERNAL_E = &PointFunctions::internal_energy;
PointProperty RHO_INTERNAL_E = &PointFunctions::rho_internal_energy;
} // namespace alias
