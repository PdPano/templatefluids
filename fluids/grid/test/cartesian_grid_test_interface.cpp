#include "cartesian_grid_test_interface.hpp"
#include "../../input_output/readers/reader.hpp"

CartesianGridTestInterface::CartesianGridTestInterface()
    : CartesianGrid::CartesianGrid()
{
}

CartesianGridTestInterface::CartesianGridTestInterface(Options& opt,
    std::istream&& initial_grid, std::istream&& mesh_details,
    std::istream&& boundary_file)
    : CartesianGrid::CartesianGrid(Reader(opt, std::move(initial_grid),
          std::move(mesh_details), std::move(boundary_file)))
{
}
