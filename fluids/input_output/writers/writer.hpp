#ifndef WRITER_HPP
#define WRITER_HPP

#include "../../grid/cartesian_grid.hpp"
#include "../../utils/point_functions.hpp"
#include "../options.hpp"
#include "default_writer.hpp"
#include "vtk_writer.hpp"
#include <iostream>
#include <string>

class Writer {
public:
    Writer(Options& opt);
    template <typename Grid>
    void write(Grid& grid, const double& t)
    {
        std::cout << std::endl
                  << "On write " << counter << "; t=" << t << std::endl;
        std::string number = std::to_string(counter);
        std::string padded_number
            = std::string(8 - number.length(), '0') + number;
        auto file_name = base_name + "_" + padded_number;

        if (output_type == "DEFAULT") {
            default_writer(grid, file_name, pf);
        }
        else if (output_type == "VTK") {
            vtk_writer(grid, file_name, pf);
        }
        else {
            std::cerr << "Output type '" << output_type << "' is not supported"
                      << std::endl;
        }
        counter++;
    }

private:
    std::string base_name;
    std::string output_type;
    int counter;
    PointFunctions pf;
};

#endif /* WRITER_HPP */
