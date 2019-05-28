#ifndef OPERATORS_OVERLOADS_HPP
#define OPERATORS_OVERLOADS_HPP
struct Point;
struct Flux;
class ShockDiscontinuity;
Point operator+(const Point& p1, const Point& p2);
Point operator-(const Point& p1, const Point& p2);
Point operator/(const Point& p1, const double& d);
Point operator*(const double& d, const Point& p1);
Point operator*(const Point& p1, const double& d);

Flux operator+(const Flux& f1, const Point& p2);
Point operator+(const Point& p1, const Flux& f2);
Flux operator*(const Flux& f1, const double& d);
Flux operator-(const Flux& f1);
Flux operator*(const double& d, const Flux& f1);
Flux operator/(const Flux& f1, const double& d);
bool operator==(const Point& lhs, const Point& rhs);
bool operator==(const ShockDiscontinuity& lhs, const ShockDiscontinuity& rhs);
#endif /* OPERATORS_OVERLOADS_HPP */
