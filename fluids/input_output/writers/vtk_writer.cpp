#include "vtk_writer.hpp"
#include "nan_checker.hpp"
#include <fstream>

void vtk_writer(CartesianGrid& grid, std::string& file_name, PointFunctions& pf)
{
    std::ofstream output(file_name + ".vtk");
    if (!output.is_open()) {
        std::cerr << "Could not open file " << file_name << std::endl;
        return;
    }

    // Header
    output << "# vtk DataFile Version 3.0" << std::endl
           << "VTK output from templateFluids" << std::endl
           << "ASCII" << std::endl
           << "DATASET RECTILINEAR_GRID" << std::endl;

    // Grid dimensions and actual coordinates
    output << "DIMENSIONS " << grid.nPointsJ << " " << grid.nPointsI << " 1"
           << std::endl;

    output << "X_COORDINATES " << grid.nPointsJ << " float" << std::endl;
    for (int j = 0; j < grid.nPointsJ; j++) {
        output << grid.X(grid.IND(0, j)) << " ";
    }
    output << std::endl;

    output << "Y_COORDINATES " << grid.nPointsI << " float" << std::endl;
    for (int i = 0; i < grid.nPointsI; i++) {
        output << grid.Y(grid.IND(i, 0)) << " ";
    }
    output << std::endl;

    output << "Z_COORDINATES " << 1 << " float" << std::endl;
    output << 0.0 << std::endl;

    // Data
    output << "POINT_DATA " << grid.nPointsTotal << std::endl;

    // Density
    output << "SCALARS Density float" << std::endl;
    output << "LOOKUP_TABLE DEFAULT" << std::endl;
    for (int i = 0; i < grid.nPointsI; i++) {
        for (int j = 0; j < grid.nPointsJ; j++) {
            int ind = grid.IND(i, j);
            output << default_if_nan(pf.rho(grid.values(ind))) << " ";
        }
        output << std::endl;
    }

    // Energy
    output << "SCALARS Energy float" << std::endl;
    output << "LOOKUP_TABLE DEFAULT" << std::endl;
    for (int i = 0; i < grid.nPointsI; i++) {
        for (int j = 0; j < grid.nPointsJ; j++) {
            int ind = grid.IND(i, j);
            output << default_if_nan(pf.e(grid.values(ind))) << " ";
        }
        output << std::endl;
    }

    // Pressure
    output << "SCALARS Pressure float" << std::endl;
    output << "LOOKUP_TABLE DEFAULT" << std::endl;
    for (int i = 0; i < grid.nPointsI; i++) {
        for (int j = 0; j < grid.nPointsJ; j++) {
            int ind = grid.IND(i, j);
            output << default_if_nan(pf.pressure(grid.values(ind))) << " ";
        }
        output << std::endl;
    }

    // Temperature
    output << "SCALARS Temperature float" << std::endl;
    output << "LOOKUP_TABLE DEFAULT" << std::endl;
    for (int i = 0; i < grid.nPointsI; i++) {
        for (int j = 0; j < grid.nPointsJ; j++) {
            int ind = grid.IND(i, j);
            output << default_if_nan(pf.temperature(grid.values(ind))) << " ";
        }
        output << std::endl;
    }

    // Mach
    output << "SCALARS Mach float" << std::endl;
    output << "LOOKUP_TABLE DEFAULT" << std::endl;
    for (int i = 0; i < grid.nPointsI; i++) {
        for (int j = 0; j < grid.nPointsJ; j++) {
            int ind = grid.IND(i, j);
            output << default_if_nan(pf.mach_number(grid.values(ind))) << " ";
        }
        output << std::endl;
    }

    // Entropy
    output << "SCALARS Entropy float" << std::endl;
    output << "LOOKUP_TABLE DEFAULT" << std::endl;
    for (int i = 0; i < grid.nPointsI; i++) {
        for (int j = 0; j < grid.nPointsJ; j++) {
            int ind = grid.IND(i, j);
            output << default_if_nan(pf.entropy(grid.values(ind))) << " ";
        }
        output << std::endl;
    }

    // Velocity
    output << "VECTORS Velocity float" << std::endl;
    for (int i = 0; i < grid.nPointsI; i++) {
        for (int j = 0; j < grid.nPointsJ; j++) {
            int ind = grid.IND(i, j);
            output << default_if_nan(pf.u(grid.values(ind))) << " "
                   << default_if_nan(pf.v(grid.values(ind))) << " 0"
                   << std::endl;
        }
    }

    if (output.fail()) {
        std::cout << "Failed while writing to file " << file_name << std::endl;
        return;
    }
}
