#ifndef GHIAS_GHOST_POINT_HPP
#define GHIAS_GHOST_POINT_HPP

struct GhiasGhostPoint {
    int ind;
    bool is_fluid_point[4];
    int neighbors_inds[4];
    double neighbors_x[4];
    double neighbors_y[4];
    double neighbors_nx[4];
    double neighbors_ny[4];
    double image_coordinate[2];
    bool is_fluid(int p) { return is_fluid_point[p]; }
};

#endif /* GHIAS_GHOST_POINT_HPP */
