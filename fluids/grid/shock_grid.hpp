/**
 * \file shock_grid.hpp
 * @brief Header for ShockGrid class
 */
#ifndef SHOCK_GRID_HPP
#define SHOCK_GRID_HPP
#include "../utils/shock_discontinuity_def.hpp"
#include "../utils/shock_discontinuity_handler.hpp"
#include "karagiozis_grid.hpp"
#include <functional>
/**
 * \class ShockGrid
 * @brief Extends KaragiozisGrid by handling moving shocks
 */
class ShockGrid : public KaragiozisGrid {
public:
    ShockGrid(Options& opt);
    ShockGrid(Reader reader, Options& opt);
    void grid_specific_pre_update(double dt) override;
    void grid_specific_update() override;
    void fill_discontinuity_map() override;
    void update_values(ShockGrid* grid_to_update_from);
    std::vector<ShockDiscontinuity>& shock_points() { return shock_points_c; }
    void compute_shock_angles();
    void extrapolate_to_shocks();
    bool merge_compatible_shocks();
    void move_shocks(double dt);
    bool update_crossed_grid_points(
        std::vector<int>& check_x_dir, std::vector<int>& check_y_dir);
    bool remove_disconnected_shocks();
    bool delete_shock_near_wall();
    void specific_print() override final;

private:
    std::vector<ShockDiscontinuity> shock_points_c;
    std::vector<ShockDiscontinuity> temp_shock_c;
    ShockHandler sh;
    bool remove_weak_shocks();
    bool detect_new_shocks(double dt);
    bool detect_new_connecting_shocks(double dt, std::vector<int>& check_x_dir,
        std::vector<int>& check_y_dir);
    void compute_shock_theta_x(ShockDiscontinuity* sp);
    void compute_shock_theta_y(ShockDiscontinuity* sp);
    void find_compatible_shocks_x(ShockDiscontinuity* sp,
        std::vector<ShockDiscontinuity>* compatible_up,
        std::vector<ShockDiscontinuity>* compatible_down);
    void find_compatible_shocks_y(ShockDiscontinuity* sp,
        std::vector<ShockDiscontinuity>* compatible_left,
        std::vector<ShockDiscontinuity>* compatible_right);
    void extend_compatible_discs(
        std::vector<ShockDiscontinuity>* compatible_discs, char compatible_side,
        int ind, DiscontinuityMap* disc_map);
    std::pair<double, double> compute_mean_position(
        const std::vector<ShockDiscontinuity>& comp_discs);
    void extrapolate_points_left(
        ShockDiscontinuity& sp, int shift_minus, DiscontinuityMap* dir_map);
    void extrapolate_points_right(ShockDiscontinuity& sp, int shift_plus,
        int shift_plus_plus, DiscontinuityMap* dir_map);
    void apply_rankine_hugoniot_jump_conditions();
    auto get_map_and_shifts(ShockDiscontinuity& sp);
    int update_crossed_by_x(ShockDiscontinuity& sp);
    int update_crossed_by_y(ShockDiscontinuity& sp);
    bool merge_shocks_in_disc_list(DiscontinuityList& disc_list);
    void merge_shocks_in_compatible_list(DiscontinuityList& compatible_list);
    void find_shocks_to_merge(DiscontinuityList* low_p_left,
        DiscontinuityList* low_p_right, DiscontinuityList& disc_list);
    auto& boundary_point_from_index(int ind);
    int update_crossed_generic_dir(ShockDiscontinuity& sp,
        const std::function<int(int)>& addInd,
        const std::function<int(int)>& subInd,
        const std::function<bool(const BoundaryPoint&)>& boundary_dir);
    bool find_wall_shock(DiscontinuityList& disc_list);
    bool create_x_shock(double dt, int ind, int ind_x);
    bool create_y_shock(double dt, int ind, int ind_y);
    const std::string base_path;
    int counter;
};
#endif /* SHOCK_GRID_HPP */
