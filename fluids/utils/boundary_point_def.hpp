#ifndef BOUNDARY_POINT_DEF_HPP
#define BOUNDARY_POINT_DEF_HPP

#include <cmath>
#include <istream>

struct BoundaryPoint {
    int ind;
    bool x_boundary, y_boundary;
    int x_type, y_type;
    double rho, u, v, T, p;
    int time_function_flag;
    double time_param1, time_param2;
    friend std::istream& operator>>(std::istream& is, BoundaryPoint& bp)
    {
        is >> bp.ind >> bp.x_boundary >> bp.y_boundary >> bp.x_type >> bp.y_type
            >> bp.rho >> bp.u >> bp.v >> bp.T >> bp.p >> bp.time_function_flag
            >> bp.time_param1 >> bp.time_param2;
        return is;
    }
    bool is_x() const { return x_boundary; }
    bool is_y() const { return y_boundary; }
    bool is_left() const { return x_type < 8 and x_boundary; }

    bool is_right() const { return x_type > 8 and x_boundary; }

    bool is_bottom() const { return y_type < 8 and y_boundary; }

    bool is_top() const { return y_type > 8 and y_boundary; }

    bool is_subsonic_inlet_x() const { return (x_type % 8) == 1; }

    bool is_subsonic_inlet_y() const { return (y_type % 8) == 1; }

    bool is_subsonic_inlet() const
    {
        return is_subsonic_inlet_x() or is_supersonic_inlet_y();
    }

    bool is_supersonic_inlet_x() const { return (x_type % 8) == 2; }

    bool is_supersonic_inlet_y() const { return (y_type % 8) == 2; }

    bool is_supersonic_inlet() const
    {
        return is_supersonic_inlet_x() or is_supersonic_inlet_y();
    }

    bool is_subsonic_outlet_x() const { return (x_type % 8) == 3; }

    bool is_subsonic_outlet_y() const { return (y_type % 8) == 3; }

    bool is_subsonic_outlet() const
    {
        return is_subsonic_outlet_x() or is_subsonic_outlet_y();
    }

    bool is_supersonic_outlet_x() const { return (x_type % 8) == 4; }

    bool is_supersonic_outlet_y() const { return (y_type % 8) == 4; }

    bool is_supersonic_outlet() const
    {
        return is_supersonic_outlet_x() or is_supersonic_outlet_y();
    }

    bool is_adiabatic_no_slip_wall_x() const { return (x_type % 8) == 5; }

    bool is_adiabatic_no_slip_wall_y() const { return (y_type % 8) == 5; }

    bool is_adiabatic_no_slip_wall() const
    {
        return is_adiabatic_no_slip_wall_x() or is_adiabatic_no_slip_wall_y();
    }

    bool is_isotermal_no_slip_wall_x() const { return (x_type % 8) == 6; }

    bool is_isotermal_no_slip_wall_y() const { return (y_type % 8) == 6; }

    bool is_isotermal_no_slip_wall() const
    {
        return is_isotermal_no_slip_wall_x() or is_isotermal_no_slip_wall_y();
    }
    double time_function(const double& t) const
    {
        switch (time_function_flag) {
        case 0:
            return 1;
        case 1:
            return linear_function(t);
        case 2:
            return exponential_function(t);
        case 3:
            return oscillatory_function(t);
        case 4:
            return multi_linear_function(t);
        default:
            throw(-1);
        };
    }

    double der_time_function(const double& t) const
    {
        const double dt = 1e-6;
        return (time_function(t + dt) - time_function(t)) / dt;
    }

    double linear_function(const double& t) const
    {
        if (t < time_param1)
            return t / time_param1;
        return 1.0;
    }

    double exponential_function(const double& t) const
    {
        return 1 - exp(-t / time_param1);
    }

    double oscillatory_function(const double& t) const
    {
        double val = sin(time_param1 * t);
        return val * val;
    }

    double multi_linear_function(const double& t) const
    {
        return linear_function(t);
    }
};

#endif /* BOUNDARY_POINT_DEF_HPP */
