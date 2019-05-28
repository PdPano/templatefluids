#ifndef CARTESIAN_GRID_TEST_INTERFACE_HPP
#define CARTESIAN_GRID_TEST_INTERFACE_HPP
#include "../cartesian_grid.hpp"

class CartesianGridTestInterface : public CartesianGrid {
public:
    CartesianGridTestInterface();
    CartesianGridTestInterface(Options& opt, std::istream&& initial_grid,
        std::istream&& mesh_details, std::istream&& boundary_file);
};

#endif /* CARTESIAN_GRID_TEST_INTERFACE_HPP */
