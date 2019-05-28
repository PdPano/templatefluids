#ifndef GRID_COMPONENTS_CONTAINER_DEF_HPP
#define GRID_COMPONENTS_CONTAINER_DEF_HPP

#include "body_discontinuity_def.hpp"
#include "boundary_point_def.hpp"
#include "ghias_ghost_point_def.hpp"
#include "point_def.hpp"
#include "shock_discontinuity_def.hpp"
#include <vector>

struct GridComponentsContainer {
    std::vector<Point> grid_c;
    std::vector<int> flags_c;
    std::vector<BoundaryPoint> boundary_c;
    std::vector<GhiasGhostPoint> ghias_points_c;
    std::vector<BodyDiscontinuity> karagiozis_points_c;
    std::vector<ShockDiscontinuity> shock_points_c;
    GridComponentsContainer()
        : ghias_points_c(std::vector<GhiasGhostPoint>(0))
        , karagiozis_points_c(std::vector<BodyDiscontinuity>(0))
        , shock_points_c(std::vector<ShockDiscontinuity>(0))
    {
    }
};

#endif /* GRID_COMPONENTS_CONTAINER_DEF_HPP */
