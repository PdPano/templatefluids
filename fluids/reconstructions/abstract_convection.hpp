#ifndef ABSTRACT_CONVECTION_HPP
#define ABSTRACT_CONVECTION_HPP

#include <memory>

class CartesianGrid;
class Derivatives;
struct PointFunctions;
struct Flux;

class Convection {
public:
    const PointFunctions& pf;
    std::shared_ptr<Derivatives> der;
    virtual Flux convection_x(const CartesianGrid& grid, int ind) const = 0;
    virtual Flux convection_y(const CartesianGrid& grid, int ind) const = 0;
    virtual void init(const CartesianGrid&) = 0;

protected:
    Convection(PointFunctions& pf_in, std::shared_ptr<Derivatives> der_in)
        : pf(pf_in)
        , der(der_in)
    {
    }
};

#endif /* ABSTRACT_CONVECTION_HPP */
