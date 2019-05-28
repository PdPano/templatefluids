#ifndef SHOCK_DISCONTINUITY_HANDLER_HPP
#define SHOCK_DISCONTINUITY_HANDLER_HPP

#include "shock_discontinuity_def.hpp"
#include <utility>

struct Point;
class Options;

class ShockHandler {
public:
    ShockHandler(Options& opt);
    ShockDiscontinuity create_shock(
        Point left_in, Point right_in, char type_in, double dx_dt = 1.0);
    void prepare_values(ShockDiscontinuity& shock);
    void compute_wall_interaction(ShockDiscontinuity& shock, double wall_theta);
    void fix_theta(ShockDiscontinuity& shock);

private:
    using Velocity = std::pair<double, double>;
    PointFunctions pf;
    const double gam;
    const double delta;
    const double angular_coeff;
    const double cfl;

    void apply_rh_from_w_mach(ShockDiscontinuity& shock, Velocity& low_vel,
        double w, double rel_mach);
    double sigma_from_mach(double mach);
    Point& low_point(ShockDiscontinuity& shock);
    Point& high_point(ShockDiscontinuity& shock);
    double estimate_rel_mach(double sigma);
    double low_c(ShockDiscontinuity& shock);
    double high_c(ShockDiscontinuity& shock);
    Velocity to_normal(Velocity vel, double theta);
    Velocity from_normal(Velocity vel, double theta);
    double u(Velocity vel);
    double v(Velocity vel);
    Velocity vel_from_point(Point& p);
    double calculate_w(double rel_mach, double low_u, double low_c);
    double calculate_rel_mach(double sigma);
    bool is_candidate(ShockDiscontinuity& shock, double dx_dt);
    double get_mach_at_wall_col(double low_u, double low_c);
};

#endif /* SHOCK_DISCONTINUITY_HANDLER_HPP */
