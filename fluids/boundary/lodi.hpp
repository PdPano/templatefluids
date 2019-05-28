/**
 * \file lodi.hpp
 * @brief Header for Lodi class
 */
#ifndef LODI_HPP
#define LODI_HPP

#include "../derivatives/derivatives.hpp"
#include <array>
#include <memory>

struct PointFunctions;
class CartesianGrid;

typedef std::array<double, 4> LodiArray;

/**
 * \class Lodi
 * @brief Class generates the Lodi relations at a grid point
 *
 * The Lodi relations are the characteristic wave amplitudes at a given
 * point and are used to define the boundary conditions
 */
class Lodi {
public:
    /**
     * @brief Constructor
     *
     * @param pf_in
     * @param der_in
     */
    Lodi(PointFunctions& pf_in, std::shared_ptr<Derivatives> der_in);

    /**
     * @brief Compute the lodi values
     *
     * @param grid Grid where the values are computed
     * @param df Derivative function (DX or DY)
     * @param vel_1 Velocity perpendicular to the boundary
     * @param vel_2 Velocity parallel to the boundary
     * @param ind Location of the boundary
     *
     * @return
     */
    LodiArray values(const CartesianGrid& grid, DerFunction df,
        alias::PointProperty vel_1, alias::PointProperty vel_2, int ind) const;

private:
    PointFunctions& pf;
    std::shared_ptr<Derivatives> der;
};

#endif /* LODI_HPP */
