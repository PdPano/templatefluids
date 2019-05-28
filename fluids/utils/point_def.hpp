#ifndef POINT_DEF_HPP
#define POINT_DEF_HPP

#include <istream>
#include <ostream>

struct Point {
    double rho_v;
    double ru_v;
    double rv_v;
    double e_v;
    Point() {}
    Point(double v1, double v2, double v3, double v4)
        : rho_v(v1)
        , ru_v(v2)
        , rv_v(v3)
        , e_v(v4)
    {
    }
    inline double rho(void) const { return rho_v; }
    inline double ru(void) const { return ru_v; }
    inline double rv(void) const { return rv_v; }
    inline double e(void) const { return e_v; }

    inline void set_rho(double val) { rho_v = val; }
    inline void set_ru(double val) { ru_v = val; }
    inline void set_rv(double val) { rv_v = val; }
    inline void set_e(double val) { e_v = val; }
    friend std::istream& operator>>(std::istream& is, Point& p)
    {
        is >> p.rho_v >> p.ru_v >> p.rv_v >> p.e_v;
        return is;
    }
    friend std::ostream& operator<<(std::ostream& os, const Point& p)
    {
        os << p.rho_v << " " << p.ru_v << " " << p.rv_v << " " << p.e_v;
        return os;
    }
};

#endif /* POINT_DEF_HPP */
