/**
 * \file irregular_derivatives.hpp
 * @brief Header for IrregularDerivatives class
 */
#ifndef IRREGULAR_DERIVATIVES_HPP
#define IRREGULAR_DERIVATIVES_HPP

class KaragiozisGrid;
class GenericDiscontinuity;
#include "regular_derivatives.hpp"

#include <functional>

typedef GenericDiscontinuity* (KaragiozisGrid::*GetDisc)(int ind) const;

/**
 * \class IrregularDerivatives
 * @brief Concrete implementation of Derivatives
 */
class IrregularDerivatives : public RegularDerivatives {

public:
    IrregularDerivatives(PointFunctions& pf_in, int shiftX_in, int shiftY_in,
        double dx_in, double dy_in);
    double DX(const CartesianGrid& grid, alias::PointProperty func,
        int ind) const override;
    double DY(const CartesianGrid& grid, alias::PointProperty func,
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

    double DXX(const CartesianGrid& grid, alias::PointProperty func,
        int ind) const override;
    double DYY(const CartesianGrid& grid, alias::PointProperty func,
        int ind) const override;
    double DXY(const CartesianGrid& grid, alias::PointProperty func,
        int ind) const override;

private:
    /**
     * @brief Holds info about discontinuities to either side
     */
    struct HasDiscs {
        bool left, right;
    };

    /**
     * @brief Returns relevant part of discontinuity point
     *
     * @param disc A discontinuity
     *
     * @return tuple containing Right point and frac distance to point at der
     */
    auto disc_point_and_frac_left(GenericDiscontinuity* disc) const;
    /**
     * @brief Returns relevant part of discontinuity point
     *
     * @param disc A discontinuity
     *
     * @return tuple containing Left point and frac distance to point at der
     */
    auto disc_point_and_frac_right(GenericDiscontinuity* disc) const;

    /**
     * @brief Checks if second derivative can be calculated
     *
     * Both parameters relate to how much room we have to take the derivative
     * @param has_discs
     * @param der_flag
     *
     * @return true if the second der is zero
     */
    bool second_der_is_zero(HasDiscs& has_discs, int der_flag) const;

    /**
     * @brief Adjust position if near the domain border
     * This will multiply the relevant shift
     *
     * @param der_flag
     *
     * @return Shift the position (+1,0,-1)
     */
    int shift_pos_to_derive(int der_flag) const;

    /**
     * @name Der_inner
     * @brief Checks if irregular or regular derivative is necessary and
     * dispatch accordingly
     *
     * @tparam Func Type of function beeing derived
     * @param grid Grid where everything is computed
     * @param func Function to derive
     * @param ind Position to derive
     *
     * @return derivative value

     * @{ */

    template <typename Func>
    auto DX_inner(const CartesianGrid& grid, Func func, int ind) const;
    template <typename Func>
    auto DY_inner(const CartesianGrid& grid, Func func, int ind) const;
    /**  @} */
    /**
     * @name Der selectors
     *
     * @tparam Func type of function beeing derived
     * @tparam first_disc Function that returns the first discontinuity to right
     * @tparam last_disc Function that return the last discontinuity to left
     * @param grid
     * @param func Function to derive
     * @param ind Position to derive
     * @param der_flag
     * @param has_discs(_close)(_far)
     * @param shift shiftX or shiftY
     * @param h dx or dy
     * @{ */

    template <typename Func, GetDisc first_disc, GetDisc last_disc>
    auto select_first_irregular_derivative(const KaragiozisGrid& grid,
        Func& func, int ind, int der_flag, HasDiscs& has_discs, int shift,
        double h) const;
    template <GetDisc first_disc, GetDisc last_disc>
    double select_second_irregular_derivative(const KaragiozisGrid& grid,
        PointToDouble& func, int ind, int der_flag, HasDiscs& has_discs_close,
        HasDiscs& has_discs_far, int shift, double h) const;

    /**  @} */

    /**
     * @name Irregular first derivatives
     * @{ */

    template <typename Func>
    auto first_der_two_discs(Func& func, GenericDiscontinuity* left,
        GenericDiscontinuity* right, double h) const;
    template <typename Func>
    auto first_der_one_disc(Func& func, const Point& q_disc, const Point& q1,
        const Point& q2, double eta, double h) const;
    /**  @} */
    /**
     * @name Irregular second derivatives
     * @{ */

    double second_der_one_disc(PointToDouble& func, const Point& q_disc,
        const Point& q1, const Point& q2, const Point& q3, double eta,
        double h) const;
    double second_der_two_disc(PointToDouble& func, const Point& q_left,
        const Point& q1, const Point& q2, const Point& q_right, double eta,
        double eps, double h) const;
    /**  @} */
    double cross_der_estimate(const CartesianGrid& grid,
        alias::PointProperty func, int ind, int i, int j) const;

    template <typename Func>
    auto zero_der() const;
};

#endif /* IRREGULAR_DERIVATIVES_HPP */
