#ifndef DERIVATIVES_HPP
#define DERIVATIVES_HPP

#include "../utils/useful_alias.hpp"
#include <functional>

class CartesianGrid;
struct Point;
struct PointFunctions;
struct Flux;

using PointToDouble = const std::function<double(const Point&)>;
using PositionToDouble = const std::function<double(int)>;

using PointToFlux = const std::function<Flux(const Point&)>;
using PositionToFlux = const std::function<Flux(int)>;

/**
 * @brief Abstract derivative class
 */
class Derivatives {
public:
    const PointFunctions& pf;
    const int shiftX;
    const int shiftY;
    const double dx;
    const double dy;

    virtual double DX(const CartesianGrid& grid, alias::PointProperty func,
        int ind) const = 0; //!< First derivative in the X direction
    virtual double DY(const CartesianGrid& grid, alias::PointProperty func,
        int ind) const = 0; //!< First derivative in the Y direction
    virtual double DXX(const CartesianGrid& grid, alias::PointProperty func,
        int ind) const = 0; //!< Second derivative in the X direction
    virtual double DYY(const CartesianGrid& grid, alias::PointProperty func,
        int ind) const = 0; //!< Second derivative in the Y direction
    virtual double DXY(const CartesianGrid& grid, alias::PointProperty func,
        int ind) const = 0; //!< Cross derivative

    virtual double DX(const CartesianGrid& grid, PointToDouble& func,
        int ind) const = 0; //!<First derivative in the X direction
    virtual double DY(const CartesianGrid& grid, PointToDouble& func,
        int ind) const = 0; //!<First derivative in the Y direction

    virtual double DX(const CartesianGrid& grid, PositionToDouble& func,
        int ind) const = 0; //!<First derivative in the X direction
    virtual double DY(const CartesianGrid& grid, PositionToDouble& func,
        int ind) const = 0; //!<First derivative in the Y direction

    virtual Flux DX(const CartesianGrid& grid, PointToFlux& func,
        int ind) const = 0; //!<First derivative in the X direction
    virtual Flux DXForward(const CartesianGrid& grid, PointToFlux& func,
        int ind) const = 0; //!< Right-sided first derivative in the X direction
    virtual Flux DXBackward(const CartesianGrid& grid, PointToFlux& func,
        int ind) const = 0; //!< Left-sided first derivative in the X direction
    virtual Flux DY(const CartesianGrid& grid, PointToFlux& func,
        int ind) const = 0; //!<First derivative in the Y direction
    virtual Flux DYForward(const CartesianGrid& grid, PointToFlux& func,
        int ind) const = 0; //!<Right-sided first derivative in the Y direction
    virtual Flux DYBackward(const CartesianGrid& grid, PointToFlux& func,
        int ind) const = 0; //!<Left-sided first derivative in the Y direction

    virtual Flux DX(const CartesianGrid& grid, PositionToFlux& func,
        int ind) const = 0; //!<First derivative in the X direction
    virtual Flux DXForward(const CartesianGrid& grid, PositionToFlux& func,
        int ind) const = 0; //!< Right-sided first derivative in the X direction
    virtual Flux DXBackward(const CartesianGrid& grid, PositionToFlux& func,
        int ind) const = 0; //!< Left-sided first derivative in the X direction
    virtual Flux DY(const CartesianGrid& grid, PositionToFlux& func,
        int ind) const = 0; //!<First derivative in the Y direction
    virtual Flux DYForward(const CartesianGrid& grid, PositionToFlux& func,
        int ind) const = 0; //!<Right-sided first derivative in the Y direction
    virtual Flux DYBackward(const CartesianGrid& grid, PositionToFlux& func,
        int ind) const = 0; //!<Left-sided first derivative in the Y direction

    virtual ~Derivatives() = default;

protected:
    /** Class is abstract, so no public constructors. Derived classes have
     * to initialize the variables*/
    Derivatives(PointFunctions& pf_in, int shiftX_in, int shiftY_in,
        double dx_in, double dy_in)
        : pf(pf_in)
        , shiftX(shiftX_in)
        , shiftY(shiftY_in)
        , dx(dx_in)
        , dy(dy_in)
    {
    }
};

typedef double (Derivatives::*DerFunction)(
    const CartesianGrid&, alias::PointProperty, int) const;

#endif /* DERIVATIVES_HPP */
