/**
 * \file ghias_shock_grid.hpp
 * @brief Header for GhiasGrid class
 */
#ifndef GHIAS_SHOCK_GRID_HPP
#define GHIAS_SHOCK_GRID_HPP

#include "ghias_grid.hpp"
#include <memory>
#include <set>
#include "../utils/shock_detectors/luisa_detector_factory.hpp"
#include "../utils/shock_detectors/luisa_shock_detector.hpp"

/**
 * \class GhiasShockGrid
 * @brief Implements the shock detector hybrid method on the ghias grid
 */
class GhiasShockGrid : public GhiasGrid {
public:
    GhiasShockGrid(Options& opt);
    GhiasShockGrid(Reader reader, Options& opt);
    virtual void grid_specific_pre_update(double) override final;
    virtual void specific_print() override final;
    const std::set<int>& to_revisit() const
    {
        return shock_detector->set_of_shocked_points();
    }

private:
    std::shared_ptr<LuisaDetector> shock_detector;
    const std::string base_path;
    int counter;
};

#endif /* GHIAS_SHOCK_GRID_HPP */
