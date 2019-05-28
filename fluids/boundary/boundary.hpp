/**
 * \file boundary.hpp
 * @brief Headers for Boundary class
 */
#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "../utils/flux_def.hpp"
#include "../utils/useful_alias.hpp"
#include "lodi.hpp"
#include <array>
#include <memory>

class CartesianGrid;
class Derivatives;
struct PointFunctions;
struct BoundaryPoint;
class DissipationTool;

/**
 * \class Boundary
 * @brief Provides methods to obtain convection, dissipation an fix the values
 * at the boundaries
 */
class Boundary {
public:
    /**
     * @brief Boundary constructor
     *
     * @param pf_in
     * @param der_in
     * @param reynolds_in
     * @param prandtl_in
     */
    Boundary(PointFunctions& pf_in, std::shared_ptr<Derivatives> der_in,
        double reynolds_in, double prandtl_in);

    /**
     * @brief Convection in x direction
     *
     * @param grid Grid where convection is computed
     * @param bp Boundary poixt where convection is computed
     * @param t Current simulation time
     *
     * @return Convection in x direction
     */
    Flux convection_x(const CartesianGrid& grid, const BoundaryPoint& bp,
        const double& t) const;
    /**
     * @brief Convection in y direction
     *
     * @param grid Grid where convection is computed
     * @param bp Boundary point where convection is computed
     * @param t Current simulation time
     *
     * @return Convection in y direction
     */
    Flux convection_y(const CartesianGrid& grid, const BoundaryPoint& bp,
        const double& t) const;

    /**
     * @brief Dissipation in x direction
     *
     * @param grid Grid where dissipation is computed
     * @param bp Boundary point where dissipation is computed
     *
     * @return Dissipation in x direction
     */
    Flux dissipation_x(
        const CartesianGrid& grid, const BoundaryPoint& bp) const;
    /**
     * @brief Dissipation in y direction
     *
     * @param grid Grid where dissipation is computed
     * @param bp Boundary point where dissipation is computed
     *
     * @return Dissipation in y direction
     */
    Flux dissipation_y(
        const CartesianGrid& grid, const BoundaryPoint& bp) const;

    /**
     * @brief Force boundary values to the ones prescribed
     *
     * @param grid Grid to fix the values
     * @param bp Boundary point to fix
     * @param t Current simulation time
     */
    void fix_boundary(
        CartesianGrid* grid, const BoundaryPoint& bp, const double& t) const;

private:
    // Utility parameters
    PointFunctions& pf;
    std::shared_ptr<Derivatives> der;
    const double gam;
    Lodi lodi;
    std::shared_ptr<DissipationTool> dissipation_tool;
    const double heat_prefactor;

    struct BoundaryConfiguration { ///< Struct for internal configuration
        const BoundaryPoint& bp;   ///< Current point
        bool is_left_or_bottom;    ///< Direction
        alias::PointProperty v1;   ///< u (if x boundary) or v (if y boundary)
        alias::PointProperty v2;   ///< v (if x boundary) or u (if y boundary)
        DerFunction df;            ///< DX (if x boundary) or DY (if y boundary)
        double t;                  ///< Current time
        BoundaryConfiguration(const BoundaryPoint& bp_in)
            : bp(bp_in)
        {
        }
    };

    typedef std::array<double, 4>
        d_array; ///< d values as defined in Poinsot & Lele

    /**
     * @brief Generates the boundary configuration structure
     *
     * @param direction 'X' or 'Y'
     * @param bp Current boundary point
     * @param t Current time
     *
     * @return Internal configuration structure
     */
    BoundaryConfiguration configure(
        char direction, const BoundaryPoint& bp, const double& t) const;

    /**
     * @name Convection functions
     * All functions change the values at lodi_a to match the boudary conditions
     * beeing applied
     * @param lodi_a Array with lodi values calculated by lodi class
     * @param grid Grid where boundary is located
     * @param bconf Local configuration struct
     * @{ */
    void subsonic_inlet(LodiArray* lodi_a, const CartesianGrid& grid,
        const BoundaryConfiguration& bconf) const;

    void supersonic_inlet(LodiArray* lodi_a, const CartesianGrid& grid,
        const BoundaryConfiguration& bconf) const;

    void subsonic_outlet(LodiArray* lodi_a, const CartesianGrid& grid,
        const BoundaryConfiguration& bconf) const;

    void supersonic_outlet(LodiArray* lodi_a, const CartesianGrid& grid,
        const BoundaryConfiguration& bconf) const;

    void adiabatic_no_slip_wall(LodiArray* lodi_a, const CartesianGrid& grid,
        const BoundaryConfiguration& bconf) const;

    void isotermal_no_slip_wall(LodiArray* lodi_a, const CartesianGrid& grid,
        const BoundaryConfiguration& bconf) const;
    /**  @} */

    /**
     * @brief Set result to correct direction
     *
     * All calculations are performed as if the boundary is in the x direction.
     * This functions resets the values to the y direction if necessary.
     *
     * @param flux_in Original result
     *
     * @return Corrected result
     */
    Flux rotate_to_y(Flux flux_in) const
    {
        return {flux_in.rho, flux_in.rv, flux_in.ru, flux_in.e};
    }

    /**
     * @brief Compute d values from lodi as in Poinsot & Lele
     *
     * @param grid Grid where the values are computed
     * @param lodi_a Lodi array computed by lodi class
     * @param bp Current boundary point
     *
     * @return d_array computed from lodi
     */
    d_array d_from_lodi(const CartesianGrid& grid, const LodiArray& lodi_a,
        const BoundaryPoint& bp) const;

    /**
     * @brief Computes the final flux from corrected lodi
     *
     * @param grid Grid where the values are computed
     * @param lodi_a Values corrected by convection functions
     * @param bconf Local configuration struct
     *
     * @return Flux computed at the boundary
     */
    Flux flux_from_lodi(const CartesianGrid& grid, const LodiArray& lodi_a,
        const BoundaryConfiguration& bconf) const;
};

#endif /* BOUNDARY_HPP */
