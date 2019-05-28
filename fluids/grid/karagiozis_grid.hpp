/**
 * \file karagiozis_grid.hpp
 * @brief Headers for KaragiozisGrid class
 */
#ifndef KARAGIOZIS_GRID_HPP
#define KARAGIOZIS_GRID_HPP

#include "../utils/body_discontinuity_def.hpp"
#include "../utils/point_functions.hpp"
#include "cartesian_grid.hpp"
#include <set>
#include <unordered_map>

class GenericDiscontinuity;

using DiscontinuityList = std::vector<GenericDiscontinuity*>;
using DiscontinuityMap = std::unordered_map<int, DiscontinuityList>;

/**
 * \class KaragiozisGrid
 * @brief Implements a grid that keeps track of body discontinuities and
 * provides methods to take special derivatives near them
 */
class KaragiozisGrid : public CartesianGrid {
public:
    KaragiozisGrid(Options& opt);
    KaragiozisGrid(Reader reader, Options& opt);

    virtual void grid_specific_update() override;

    /**
     * @name Accessors
     * @{ */

    GenericDiscontinuity* first_disc_x(int ind) const;
    GenericDiscontinuity* last_disc_x(int ind) const;
    GenericDiscontinuity* first_disc_y(int ind) const;
    GenericDiscontinuity* last_disc_y(int ind) const;

    DiscontinuityMap* discontinuity_map_x() { return &disc_map_x; }
    DiscontinuityMap* discontinuity_map_y() { return &disc_map_y; }

    double disc_X(const GenericDiscontinuity* disc)
    {
        if (disc->is_x()) {
            return CartesianGrid::X(disc->ind) + dx * disc->frac;
        }
        return CartesianGrid::X(disc->ind);
    }
    double disc_Y(const GenericDiscontinuity* disc)
    {
        if (disc->is_y()) {
            return CartesianGrid::Y(disc->ind) + dy * disc->frac;
        }
        return CartesianGrid::Y(disc->ind);
    }
    auto& body_points() { return karagiozis_points_c; }

    /**  @} */

    bool has_discont_x(int ind, int shift = 0) const;
    bool has_discont_y(int ind, int shift = 0) const;
    void clear_discontinuity_map();
    virtual void fill_discontinuity_map();

    const std::set<int>& to_revisit() const { return points_to_revisit; }
    PointFunctions pf;

protected:
    void sort_lists(); ///< Sorts the lists inside DiscontinuityMap
    void set_points_to_revisit(GenericDiscontinuity* disc);
    void safe_insert_to_points_to_revisit(int ind);
    std::set<int> points_to_revisit; ///< Points close to discontinuities

private:
    std::vector<BodyDiscontinuity> karagiozis_points_c;
    DiscontinuityMap disc_map_x;
    DiscontinuityMap disc_map_y;
    const DiscontinuityMap* discontinuity_map_x() const { return &disc_map_x; }
    const DiscontinuityMap* discontinuity_map_y() const { return &disc_map_y; }

    /**
     * @name Extrapolation functions
     * @{ */
    /**
     * Returns the relevant DiscontinuityMap and the three shifts needed to
     * compute the extrapolation from the left and from the right of the
     * discontinuity
     */
    auto get_map_and_shifts(BodyDiscontinuity& bd);
    void extrapolate_left(
        BodyDiscontinuity& bd, int shifted_ind, DiscontinuityMap* dir_map);
    void extrapolate_right(BodyDiscontinuity& bd, int shift_plus,
        int shift_plus_plus, DiscontinuityMap* dir_map);

    double linear_extrapolation(double a, double b, double eps)
    {
        // with limiter
        // return b + (b - a) * eps * (1 - eps);
        // no limiter
        return b + (b - a) * eps;
    }

    /**  @} */
};

#endif /* KARAGIOZIS_GRID_HPP */
