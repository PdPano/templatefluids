/**
 * \file cartesian_grid.cpp
 * @brief Implementation of CartesianGrid class
 */
#include "./cartesian_grid.hpp"
#include "../input_output/options.hpp"
#include "../input_output/readers/reader.hpp"
#include "../utils/flag_handler.hpp"

CartesianGrid::CartesianGrid()
    : nPointsI(10)
    , nPointsJ(10)
    , nPointsTotal(nPointsI * nPointsJ)
    , dx(0.1)
    , dy(0.1)
    , xmin(0)
    , ymin(0)
    , points_c(std::vector<Point>(nPointsTotal))
{
}

CartesianGrid::CartesianGrid(Options& opt)
    : CartesianGrid(Reader(opt))
{
}

CartesianGrid::CartesianGrid(Reader reader)
    : nPointsI(reader.nPointsI())
    , nPointsJ(reader.nPointsJ())
    , nPointsTotal(reader.nPointsTotal())
    , dx(reader.dx())
    , dy(reader.dy())
    , xmin(reader.xmin())
    , ymin(reader.ymin())
    , points_c(reader.grid())
    , flags_c(reader.flags())
    , boundary_c(reader.boundary())
{
}

void CartesianGrid::update_values(CartesianGrid* grid_to_update_from)
{
    bool should_throw = false;
    if (grid_to_update_from->nPointsI != nPointsI
        or grid_to_update_from->nPointsJ != nPointsJ) {
        std::cerr << "Incompatible dimensions during copy!!!" << std::endl;
        should_throw = true;
    }
    if (grid_to_update_from->dx != dx or grid_to_update_from->dy != dy) {
        std::cerr << "Incompatible grid spacing during copy!!!" << std::endl;
        should_throw = true;
    }
    if (grid_to_update_from->xmin != xmin
        or grid_to_update_from->ymin != ymin) {
        std::cerr << "Incompatible xmin/ymin during copy!!!" << std::endl;
        should_throw = true;
    }
    if (should_throw) {
        throw(-1);
    }
    points_c = grid_to_update_from->points_c;
}

bool CartesianGrid::ind_is_valid(int ind) const
{
    return (ind >= 0 and ind < nPointsI * nPointsJ
        and flag_functions::point_type(flags_c[ind]) != SOLID_POINT);
}
