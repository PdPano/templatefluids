/**
 * \file cartesian_grid.hpp
 * @brief Header for CartesianGrid class
 */
#ifndef CARTESIAN_GRID_HPP
#define CARTESIAN_GRID_HPP

class Options;
class Reader;
#include "../utils/boundary_point_def.hpp"
#include "../utils/point_def.hpp"

#include <utility>
#include <vector>

/**
 * \class CartesianGrid
 * @brief Basic grid providing info common to all grids
 *
 * This is by far the most important grid in the code base. Every calculation is
 * performed on a grid, and all grids derive from this one. Be careful changing
 * things here and, if you do make changes, always run the tests to make sure
 * nothing broke elsewhere.
 * The grid is a rectangular Cartesian mesh with fixed spacing in both
 * direcitons
 */
class CartesianGrid {
public:
    CartesianGrid(Options& opt);                        ///< Main constructor
    CartesianGrid(const CartesianGrid& grid) = default; ///< Copy constructor
    const int nPointsI;     ///< Number of grid lines
    const int nPointsJ;     ///< Number of grid columns
    const int nPointsTotal; ///< nPointsI*nPointsJ
    const double dx;        ///< Space between columns
    const double dy;        ///< Space between lines
    const double xmin;      ///< Bottom left x coordinate
    const double ymin;      ///< Bottom left y coordinate

    /**
     */

    /**
     * @brief Copies all point values to this grid
     *
     * @param grid_to_update_from
     */
    void update_values(CartesianGrid* grid_to_update_from);

    /**
     * @name Accessors
     * @{ */
    inline const Point& values(int ind) const { return points_c[ind]; }
    inline double rho(int ind) const { return values(ind).rho(); }
    inline double ru(int ind) const { return values(ind).ru(); }
    inline double rv(int ind) const { return values(ind).rv(); }
    inline double e(int ind) const { return values(ind).e(); }
    inline int flag(int ind) const { return flags_c[ind]; }

    inline const std::vector<BoundaryPoint>& boundary(void) const
    {
        return boundary_c;
    }

    double X(int ind) const { return xmin + dx * indJ(ind); }
    double Y(int ind) const { return ymin + dy * indI(ind); }
    /**  @} */

    /**
     * @name Setters
     * @{ */

    inline void setRho(double val, int ind) { _values(ind).set_rho(val); }
    inline void setRU(double val, int ind) { _values(ind).set_ru(val); }
    inline void setRV(double val, int ind) { _values(ind).set_rv(val); }
    inline void setE(double val, int ind) { _values(ind).set_e(val); }
    inline void set_values(Point p, int ind) { _values(ind) = p; }
    /**  @} */

    /**
     * @name Updates during time integration
     * @{ */

    /**
     * @brief Executed every RK level
     */
    virtual void grid_specific_update() {}

    /**
     * @brief Executed once before a RK step
     *
     * @param dt Current time step
     */
    virtual void grid_specific_pre_update(double /*dt*/) {}

    /**
     * @brief Executed once after a RK step
     *
     * @param dt Current time step
     */
    virtual void grid_specific_pos_update(double /*dt*/) {}
    /**  @} */

    /**
     * @name Index manipulation functions
     * @{ */

    int indI(int ind) const { return ind / nPointsJ; /*C-style arrays*/ }
    int indJ(int ind) const { return ind % nPointsJ; /*C-style arrays*/ }
    int indIPlusOne(int ind) const { return IND(indI(ind) + 1, indJ(ind)); }
    int indIMinusOne(int ind) const { return IND(indI(ind) - 1, indJ(ind)); }
    int indJPlusOne(int ind) const { return IND(indI(ind), indJ(ind) + 1); }
    int indJMinusOne(int ind) const { return IND(indI(ind), indJ(ind) - 1); }
    bool i_is_valid(int i) const { return (i >= 0 and i < nPointsI); }
    bool j_is_valid(int j) const { return (j >= 0 and j < nPointsJ); }
    bool ind_is_valid(int ind) const;
    int IND(int i, int j) const //!< Boundary-aware index calculation
    {
        if (i_is_valid(i) and j_is_valid(j))
            return i * nPointsJ + j;
        else
            return -1;
    }
    int IND(std::pair<int, int> p) const { return IND(p.first, p.second); }

    int shiftX(void) const { return 1; /*C-style arrays*/ }
    int shiftY(void) const { return nPointsJ; /*C-style arrays*/ }

    std::pair<int, int> indIJ(int ind) const
    {
        return std::make_pair(indI(ind), indJ(ind));
    }
    /**  @} */
    virtual void specific_print() {}

protected:
    CartesianGrid();
    CartesianGrid(Reader reader);

private:
    std::vector<Point> points_c;
    std::vector<int> flags_c;
    std::vector<BoundaryPoint> boundary_c;
    inline Point& _values(int ind) { return points_c[ind]; }
};

#endif
