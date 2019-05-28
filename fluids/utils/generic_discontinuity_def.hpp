#ifndef GENERIC_DISCONTINUITY_DEF_HPP
#define GENERIC_DISCONTINUITY_DEF_HPP

#include "point_def.hpp"
class GenericDiscontinuity {
public:
    GenericDiscontinuity(Point left_in, Point right_in)
        : theta(0)
        , left_p(left_in)
        , right_p(right_in)
    {
    }
    virtual ~GenericDiscontinuity() = default;
    char type;
    int ind;
    double frac;
    double theta;
    bool is_x() const { return (type == 'x'); }
    bool is_y() const { return (type == 'y'); }
    Point& left() { return left_p; }
    Point& right() { return right_p; }
    const Point& cleft() const { return left_p; }
    const Point& cright() const { return right_p; }

private:
    Point left_p;
    Point right_p;
};

#endif /* GENERIC_DISCONTINUITY_DEF_HPP */
