#ifndef FLUX_DEF_HPP
#define FLUX_DEF_HPP

#include <ostream>

struct Flux {
    double rho;
    double ru;
    double rv;
    double e;
    Flux operator+(const Flux& rhs)
    {
        return {rho + rhs.rho, ru + rhs.ru, rv + rhs.rv, e + rhs.e};
    }
    Flux operator-(const Flux& rhs)
    {
        return {rho - rhs.rho, ru - rhs.ru, rv - rhs.rv, e - rhs.e};
    }

    Flux& operator+=(const Flux& rhs)
    {
        rho += rhs.rho;
        ru += rhs.ru;
        rv += rhs.rv;
        e += rhs.e;
        return *this;
    }
    friend std::ostream& operator<<(std::ostream& os, const Flux& p)
    {
        os << p.rho << " " << p.ru << " " << p.rv << " " << p.e;
        return os;
    }
};

#endif /* FLUX_DEF_HPP */
