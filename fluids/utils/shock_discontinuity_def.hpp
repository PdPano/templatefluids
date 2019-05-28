#ifndef SHOCK_DISCONTINUITY_DEF_HPP
#define SHOCK_DISCONTINUITY_DEF_HPP

#include "generic_discontinuity_def.hpp"
#include "point_functions.hpp"
#include <cmath>
#include <ostream>
class ShockDiscontinuity : public GenericDiscontinuity {
public:
    ShockDiscontinuity(Point left_in, Point right_in, char type_in)
        : GenericDiscontinuity(left_in, right_in)
        , is_connected(true)
    {
        type = type_in;
        if (type == 'x') {
            theta = 0.0;
        }
        else {
            // if (type == 'y') {
            theta = M_PI / 2.;
        }
    }
    ShockDiscontinuity()
        : ShockDiscontinuity({1.0, 0, 0, 2.0}, {1.0, 0, 0, 2.0}, 'x')
    {
        sigma = 1.0;
    }

    char low_pressure_side;
    double w;
    double sigma;
    bool is_connected;

    bool is_strong() { return sigma >= 1.05; }
    bool is_weak() { return !(is_strong()); }
    friend std::ostream& operator<<(std::ostream& os, ShockDiscontinuity& sd)
    {
        auto& left_p = sd.left();
        auto& right_p = sd.right();
        os << "left: " << left_p.rho() << " " << left_p.ru() << " "
           << left_p.rv() << " " << left_p.e() << std::endl;
        os << "right: " << right_p.rho() << " " << right_p.ru() << " "
           << right_p.rv() << " " << right_p.e() << std::endl;
        os << "ind: " << sd.ind << " frac: " << sd.frac << std::endl;
        os << "w: " << sd.w << " low_pressure_side: " << sd.low_pressure_side
           << std::endl;
        os << "sigma: " << sd.sigma << " is_weak(): " << sd.is_weak()
           << std::endl;
        os << "type: " << sd.type << std::endl;
        os << "theta: " << sd.theta << std::endl;
        os << "is_connected: " << sd.is_connected << std::endl;
        return os;
    }
};

#endif /* SHOCK_DISCONTINUITY_DEF_HPP */
