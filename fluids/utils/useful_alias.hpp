#ifndef USEFUL_ALIAS_HPP
#define USEFUL_ALIAS_HPP

struct Point;
struct PointFunctions;
namespace alias {

typedef double (PointFunctions::*PointProperty)(const Point& p) const;

extern PointProperty RHO;
extern PointProperty RU;
extern PointProperty RV;
extern PointProperty E;
extern PointProperty P;
extern PointProperty C;
extern PointProperty T;
extern PointProperty RUV;
extern PointProperty U;
extern PointProperty V;
extern PointProperty RU2;
extern PointProperty RV2;
extern PointProperty INTERNAL_E;
extern PointProperty RHO_INTERNAL_E;
}

#endif /* USEFUL_ALIAS_HPP */
