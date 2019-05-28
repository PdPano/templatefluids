/**
 * \file ghias_shock_grid.cpp
 * @brief Implementation of GhiasShockGrid class
 */
#include "ghias_shock_grid.hpp"
#include "../input_output/options.hpp"
#include "../input_output/readers/reader.hpp"

GhiasShockGrid::GhiasShockGrid(Options& opt)
    : GhiasShockGrid(Reader(opt), opt)
{
}

GhiasShockGrid::GhiasShockGrid(Reader reader, Options& opt)
    : GhiasGrid(reader, opt)
    , shock_detector(
          create_luisa_detector(opt, reader.nPointsI(), reader.nPointsJ()))
    , base_path(opt.output_base_path())
    , counter(opt.output_counter())
{
}

void GhiasShockGrid::grid_specific_pre_update(double /*unused*/)
{
    shock_detector->detect_shocks(*this);
}

void GhiasShockGrid::specific_print()
{
    std::string number = std::to_string(counter);
    std::string padded_number = std::string(8 - number.length(), '0') + number;
    auto file_name = base_path + "discontinuity" + "_" + padded_number + ".txt";
    std::ofstream output(file_name);
    output << "x,y" << std::endl;
    for (auto& ind : to_revisit()) {
        output << X(ind) << "," << Y(ind) << std::endl;
    }
    counter++;
}
