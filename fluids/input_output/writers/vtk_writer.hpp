#ifndef VTK_WRITER_HPP
#define VTK_WRITER_HPP

#include "../../grid/cartesian_grid.hpp"
#include "../../utils/point_functions.hpp"
#include <iomanip>
#include <iostream>
#include <string>

void vtk_writer(
    CartesianGrid& grid, std::string& file_name, PointFunctions& pf);

#endif /* VTK_WRITER_HPP */
