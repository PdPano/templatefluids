#ifndef READER_HPP
#define READER_HPP

#include "default_reader.hpp"

#include "../../utils/body_discontinuity_def.hpp"
#include "../../utils/boundary_point_def.hpp"
#include "../../utils/ghias_ghost_point_def.hpp"
#include "../../utils/grid_components_container_def.hpp"
#include "../../utils/grid_constants_container_def.hpp"
#include "../../utils/point_def.hpp"
#include "../options.hpp"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

class Reader {
public:
    Reader(Options& opt);
    Reader(Options& opt, std::istream&& initial_conditions,
        std::istream&& mesh_details, std::istream&& boundary_file,
        std::istream&& immersed_interface = std::move(std::istringstream("")),
        std::istream&& shock_file = std::move(std::istringstream("")));
    std::vector<Point> grid(void) { return local_components.grid_c; }
    std::vector<int> flags(void) { return local_components.flags_c; }
    std::vector<BoundaryPoint> boundary(void)
    {
        return local_components.boundary_c;
    }
    std::vector<GhiasGhostPoint> ghias_ghost_points()
    {
        return local_components.ghias_points_c;
    }
    std::vector<BodyDiscontinuity> karagiozis_body_points()
    {
        return local_components.karagiozis_points_c;
    }
    std::vector<ShockDiscontinuity> shock_points()
    {
        return local_components.shock_points_c;
    }
    int nPointsI(void) { return local_container.nPointsI; }
    int nPointsJ(void) { return local_container.nPointsJ; }
    int nPointsTotal(void) { return local_container.nPointsTotal; }
    double dx(void) { return local_container.dx; }
    double dy(void) { return local_container.dy; }
    double xmin(void) { return local_container.xmin; }
    double ymin(void) { return local_container.ymin; }

private:
    GridComponentsContainer local_components;
    GridConstantsContainer local_container;
};

#endif /* READER_HPP */
