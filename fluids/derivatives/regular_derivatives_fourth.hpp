/**
 * \file regular_derivatives_fourth.hpp
 * @brief Header for RegularDerivativesFourth class
 */
#ifndef REGULAR_DERIVATIVES_FOURTH_HPP
#define REGULAR_DERIVATIVES_FOURTH_HPP

#include "derivatives.hpp"

/**
 * \class RegularDerivativesFourth
 * @brief Concrete implementation of Derivatives
 */
class RegularDerivativesFourth : public Derivatives {

public:
    RegularDerivativesFourth(PointFunctions& pf_in, int shiftX_in,
        int shiftY_in, double dx_in, double dy_in);

    double DX(const CartesianGrid& grid, alias::PointProperty func,
        int ind) const override;
    double DY(const CartesianGrid& grid, alias::PointProperty func,
        int ind) const override;
    double DXX(const CartesianGrid& grid, alias::PointProperty func,
        int ind) const override;
    double DYY(const CartesianGrid& grid, alias::PointProperty func,
        int ind) const override;
    double DXY(const CartesianGrid& grid, alias::PointProperty func,
        int ind) const override;
    double DX(
        const CartesianGrid& grid, PointToDouble& func, int ind) const override;
    double DY(
        const CartesianGrid& grid, PointToDouble& func, int ind) const override;
    double DX(const CartesianGrid& grid, PositionToDouble& func,
        int ind) const override;
    double DY(const CartesianGrid& grid, PositionToDouble& func,
        int ind) const override;
    Flux DX(
        const CartesianGrid& grid, PointToFlux& func, int ind) const override;
    Flux DXForward(
        const CartesianGrid& grid, PointToFlux& func, int ind) const override;
    Flux DXBackward(
        const CartesianGrid& grid, PointToFlux& func, int ind) const override;
    Flux DY(
        const CartesianGrid& grid, PointToFlux& func, int ind) const override;
    Flux DYForward(
        const CartesianGrid& grid, PointToFlux& func, int ind) const override;
    Flux DYBackward(
        const CartesianGrid& grid, PointToFlux& func, int ind) const override;

    Flux DX(const CartesianGrid& grid, PositionToFlux& func,
        int ind) const override; //!<First derivative in the X direction
    Flux DXForward(const CartesianGrid& grid, PositionToFlux& func,
        int ind)
        const override; //!< Right-sided first derivative in the X direction
    Flux DXBackward(const CartesianGrid& grid, PositionToFlux& func,
        int ind)
        const override; //!< Left-sided first derivative in the X direction
    Flux DY(const CartesianGrid& grid, PositionToFlux& func,
        int ind) const override; //!<First derivative in the Y direction
    Flux DYForward(const CartesianGrid& grid, PositionToFlux& func,
        int ind)
        const override; //!<Right-sided first derivative in the Y direction
    Flux DYBackward(const CartesianGrid& grid, PositionToFlux& func,
        int ind)
        const override; //!<Left-sided first derivative in the Y direction

private:
    /**
     * @name Derivatives selector
     * @{ */

    /**
     * @brief First derivative dispatcher
     *
     * @tparam Func Type of function beeing derived
     * @param func Ref to function beeing derived
     * @param ind Where to compute
     * @param shift How big is the jump between indexes (different in x and y)
     * @param h grid spacing (dx or dy)
     * @param der_flag How much grid space to either side
     *
     * @return First derivative value (double or flux)
     */
    template <typename Func>
    auto first_der_selector(
        Func& func, int ind, int shift, double h, uint16_t der_flag) const;
    /**
     * @brief Second derivative dispatcher
     *
     * @param func Ref to function beeing derived
     * @param ind Where to compute
     * @param shift How big is the jump between indexes (different in x and y)
     * @param h grid spacing (dx or dy)
     * @param der_flag How much grid space to either side
     *
     * @return Second derivative value (double)
     */
    double second_der_selector(PositionToDouble& func, int ind, int shift,
        double h, uint16_t der_flag) const;

    /**  @} */
    /**
     * @name First derivative types
     * Types of regular derivatives. Accessed via first_der_selector()
     * @{ */

    template <typename Func>
    auto first_derivative_central(
        Func& func, int ind, int shift, double h) const;
    template <typename Func>
    auto first_derivative_right(Func& func, int ind, int shift, double h) const;
    template <typename Func>
    auto first_derivative_left(Func& func, int ind, int shift, double h) const;
    template <typename Func>
    auto first_derivative_almost_left(
        Func& func, int ind, int shift, double h) const;
    template <typename Func>
    auto first_derivative_almost_right(
        Func& func, int ind, int shift, double h) const;

    /**  @} */
    /**
     * @name Second derivative types
     * Types of regular second derivatives. Accesed via second_der_selector()
     * @{ */

    double second_derivative_central(
        PositionToDouble& func, int ind, int shift, double h) const;
    double second_derivative_right(
        PositionToDouble& func, int ind, int shift, double h) const;
    double second_derivative_left(
        PositionToDouble& func, int ind, int shift, double h) const;
    double second_derivative_almost_right(
        PositionToDouble& func, int ind, int shift, double h) const;
    double second_derivative_almost_left(
        PositionToDouble& func, int ind, int shift, double h) const;
    /**  @} */
    template <typename Func>
    auto zero_der() const;
};

#endif /* REGULAR_DERIVATIVES_FOURTH_HPP */
