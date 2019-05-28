#include "default_writer.hpp"
#include "../../grid/cartesian_grid.hpp"
#include "../../grid/karagiozis_grid.hpp"
#include "../../utils/point_functions.hpp"
#include "nan_checker.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

void default_writer(
    CartesianGrid& grid, std::string& file_name, PointFunctions& pf)
{
    std::ofstream output(file_name + ".txt");

    output << "x,y,rho,u,v,E,T,p,mach,entropy" << std::endl;
    for (int ind = 0; ind < grid.nPointsTotal; ind++) {
        auto p = grid.values(ind);
        output << std::scientific << std::setprecision(15) << grid.X(ind) << ","
               << grid.Y(ind) << "," << default_if_nan(pf.rho(p)) << ","
               << default_if_nan(pf.u(p)) << "," << default_if_nan(pf.v(p))
               << "," << default_if_nan(pf.e(p)) << ","
               << default_if_nan(pf.temperature(p)) << ","
               << default_if_nan(pf.pressure(p)) << ","
               << default_if_nan(pf.mach_number(p)) << ","
               << default_if_nan(pf.entropy(p)) << std::endl;
    }
    /*    try {
            auto& local_grid = dynamic_cast<KaragiozisGrid&>(grid);
            default_writer(local_grid, output, pf);
        }
        catch (...) {
        }*/
}

void default_writer(
    KaragiozisGrid& grid, std::ofstream& output, PointFunctions& pf)
{
    for (auto& disc : grid.body_points()) {
        auto p = disc.left();
        output << std::scientific << std::setprecision(15) << grid.disc_X(&disc)
               << "," << grid.disc_Y(&disc) << "," << pf.rho(p) << ","
               << pf.u(p) << "," << pf.v(p) << "," << pf.e(p) << ","
               << pf.temperature(p) << "," << pf.pressure(p) << ","
               << pf.mach_number(p) << "," << pf.entropy(p) << std::endl;
        p = disc.right();
        output << std::scientific << std::setprecision(15) << grid.disc_X(&disc)
               << "," << grid.disc_Y(&disc) << "," << pf.rho(p) << ","
               << pf.u(p) << "," << pf.v(p) << "," << pf.e(p) << ","
               << pf.temperature(p) << "," << pf.pressure(p) << ","
               << pf.mach_number(p) << "," << pf.entropy(p) << std::endl;
    }
}
