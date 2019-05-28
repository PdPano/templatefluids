#ifndef POINT_FUNCTIONS_HPP
#define POINT_FUNCTIONS_HPP

struct Point;

struct PointFunctions {
    PointFunctions(double mach_in, double gam_in);
    const double mach, gam;

    double rho(const Point& p) const;
    double ru(const Point& p) const;
    double rv(const Point& p) const;
    double e(const Point& p) const;
    double pressure(const Point& p) const;
    double sound_speed(const Point& p) const;
    double temperature(const Point& p) const;
    double ruv(const Point& p) const;
    double u(const Point& p) const;
    double v(const Point& p) const;
    double ru2(const Point& p) const;
    double rv2(const Point& p) const;
    double internal_energy(const Point& p) const;
    double rho_internal_energy(const Point& p) const;
    double e_from_T(const Point& p, const double& T) const;
    double entropy(const Point& p) const;
    double mach_number(const Point& p) const;
};

#endif /* POINT_FUNCTIONS_HPP */
