/**
 * \file shock_grid.cpp
 * @brief Implementation of ShockGrid class
 */
#include "shock_grid.hpp"
#include "../input_output/readers/reader.hpp"
#include "../utils/debug_def.hpp"
#include "../utils/flag_handler.hpp"
#include "../utils/operators_overloads.hpp"
#include <algorithm>
#include <iostream>
ShockGrid::ShockGrid(Options& opt)
    : ShockGrid(Reader(opt), opt)
{
}

ShockGrid::ShockGrid(Reader reader, Options& opt)
    : KaragiozisGrid(reader, opt)
    , shock_points_c(reader.shock_points())
    , sh(ShockHandler(opt))
    , base_path(opt.output_base_path())
    , counter(opt.output_counter())
{
    ShockGrid::fill_discontinuity_map();
    std::ofstream outfile("./output/shocks.txt");
    for (auto& sp : shock_points_c) {
        outfile << disc_X(&sp) << "," << disc_Y(&sp) << std::endl;
    }
}

auto ShockGrid::get_map_and_shifts(ShockDiscontinuity& sp)
{
    int shift_minus;
    int shift_plus;
    int shift_plus_plus;
    if (sp.is_x()) {
        shift_minus = indJMinusOne(sp.ind);
        shift_plus = indJPlusOne(sp.ind);
        shift_plus_plus = indJPlusOne(shift_plus);
        if (!ind_is_valid(shift_minus)) {
            shift_minus = sp.ind;
        }
        if (!ind_is_valid(shift_plus_plus)) {
            shift_plus_plus = shift_plus;
        }
        return std::make_tuple(
            discontinuity_map_x(), shift_minus, shift_plus, shift_plus_plus);
    }
    // Assuming sp is y
    shift_minus = indIMinusOne(sp.ind);
    shift_plus = indIPlusOne(sp.ind);
    shift_plus_plus = indIPlusOne(shift_plus);
    if (!ind_is_valid(shift_minus)) {
        shift_minus = sp.ind;
    }
    if (!ind_is_valid(shift_plus_plus)) {
        shift_plus_plus = shift_plus;
    }
    return std::make_tuple(
        discontinuity_map_y(), shift_minus, shift_plus, shift_plus_plus);
}

void ShockGrid::fill_discontinuity_map()
{
    points_to_revisit.clear();
    for (auto& sd : shock_points_c) {
        DiscontinuityMap* dir_map;
        if (sd.is_x()) {
            dir_map = discontinuity_map_x();
        }
        else {
            dir_map = discontinuity_map_y();
        }
        (*dir_map)[sd.ind].push_back(&sd);
        set_points_to_revisit(&sd);
    }
    KaragiozisGrid::fill_discontinuity_map();
    sort_lists();
}

void ShockGrid::grid_specific_update()
{
    KaragiozisGrid::grid_specific_update();
    extrapolate_to_shocks();
    apply_rankine_hugoniot_jump_conditions();
}

void ShockGrid::grid_specific_pre_update(double dt)
{
    std::vector<int> check_x_dir;
    std::vector<int> check_y_dir;
    // Detect shocks
    bool shocks_detected = detect_new_shocks(dt);
    if (shocks_detected) {
        clear_discontinuity_map();
        ShockGrid::fill_discontinuity_map();
    }

    compute_shock_angles();

    extrapolate_to_shocks();

    apply_rankine_hugoniot_jump_conditions();

    move_shocks(dt);

    bool crossed_point = update_crossed_grid_points(check_x_dir, check_y_dir);
    bool shocks_merged = merge_compatible_shocks();
    bool weak_shocks_removed = remove_weak_shocks();
    bool shock_near_wall_deleted = delete_shock_near_wall();
    bool disconnected_shocks_removed = remove_disconnected_shocks();

    if (crossed_point or shocks_merged or weak_shocks_removed
        or disconnected_shocks_removed or shock_near_wall_deleted) {
        clear_discontinuity_map();
        ShockGrid::fill_discontinuity_map();
    }

    bool detected_connecting
        = detect_new_connecting_shocks(dt, check_x_dir, check_y_dir);
    if (detected_connecting) {
        clear_discontinuity_map();
        ShockGrid::fill_discontinuity_map();
    }
}

auto& ShockGrid::boundary_point_from_index(int ind)
{
    for (auto& bp : boundary()) {
        if (bp.ind == ind) {
            return bp;
        }
    }
    throw(-1);
}

bool ShockGrid::update_crossed_grid_points(
    std::vector<int>& check_x_dir, std::vector<int>& check_y_dir)
{
    // This is not ready yet
    // Must compute weighted average when two discontinuities crosses the same
    // point
    // [DONE] and also check if the point crossed is an outlet!
    int crossed_index;
    for (auto& sp : shock_points_c) {
        if (sp.is_x()) {
            crossed_index = update_crossed_by_x(sp);
            if (crossed_index >= 0) {
                check_y_dir.push_back(crossed_index);
            }
        }
        else {
            crossed_index = update_crossed_by_y(sp);
            if (crossed_index >= 0) {
                check_x_dir.push_back(crossed_index);
            }
        }
    }

    if (check_x_dir.size() + check_y_dir.size() > 0u) {
        return true;
    }
    return false;
}

int ShockGrid::update_crossed_by_x(ShockDiscontinuity& sp)
{
    auto addInd = [&](int ind) { return indJPlusOne(ind); };
    auto subInd = [&](int ind) { return indJMinusOne(ind); };
    auto boundary_dir = [](const BoundaryPoint& bp) { return bp.is_x(); };
    return update_crossed_generic_dir(sp, addInd, subInd, boundary_dir);
}

int ShockGrid::update_crossed_by_y(ShockDiscontinuity& sp)
{
    auto addInd = [&](int ind) { return indIPlusOne(ind); };
    auto subInd = [&](int ind) { return indIMinusOne(ind); };
    auto boundary_dir = [](const BoundaryPoint& bp) { return bp.is_y(); };
    return update_crossed_generic_dir(sp, addInd, subInd, boundary_dir);
}

int ShockGrid::update_crossed_generic_dir(ShockDiscontinuity& sp,
    const std::function<int(int)>& addInd,
    const std::function<int(int)>& subInd,
    const std::function<bool(const BoundaryPoint&)>& boundary_dir)
{
    int crossed_something = -1;
    if (sp.frac < 0) {
        crossed_something = sp.ind;
        sp.frac += 1;
        set_values(sp.right(), sp.ind);
        if (flag_functions::point_type(flag(sp.ind)) == BOUNDARY_POINT) {
            auto& bp = boundary_point_from_index(sp.ind);
            if (boundary_dir(bp)) {
                sp.is_connected = false;
            }
        }
        sp.ind = subInd(sp.ind);
    }
    else if (sp.frac >= 1) {
        sp.frac -= 1;
        sp.ind = addInd(sp.ind);
        set_values(sp.left(), sp.ind);
        if (flag_functions::point_type(flag(sp.ind)) == BOUNDARY_POINT) {
            auto& bp = boundary_point_from_index(sp.ind);
            if (boundary_dir(bp)) {
                sp.is_connected = false;
            }
        }
        crossed_something = sp.ind;
    }
    if (sp.frac < 0 or sp.frac >= 1) {
        sp.is_connected = false;
    }
    return crossed_something;
}

bool ShockGrid::merge_compatible_shocks()
{
    bool merged_shocks = false;
    auto* disc_map = discontinuity_map_x();
    for (auto& disc_list : (*disc_map)) {
        merged_shocks |= merge_shocks_in_disc_list(disc_list.second);
    }

    disc_map = discontinuity_map_y();
    for (auto& disc_list : (*disc_map)) {
        merged_shocks |= merge_shocks_in_disc_list(disc_list.second);
    }

    shock_points_c.insert(
        shock_points_c.end(), temp_shock_c.begin(), temp_shock_c.end());
    temp_shock_c.clear();

    return merged_shocks;
}

bool ShockGrid::merge_shocks_in_disc_list(DiscontinuityList& disc_list)
{
    bool merged_shocks = false;
    DiscontinuityList low_pressure_side_left(0);
    DiscontinuityList low_pressure_side_right(0);
    find_shocks_to_merge(
        &low_pressure_side_left, &low_pressure_side_right, disc_list);

    if (low_pressure_side_left.size() > 1u) {
        merge_shocks_in_compatible_list(low_pressure_side_left);
        merged_shocks = true;
    }

    if (low_pressure_side_right.size() > 1u) {
        merge_shocks_in_compatible_list(low_pressure_side_right);
        merged_shocks = true;
    }
    return merged_shocks;
}

void ShockGrid::merge_shocks_in_compatible_list(
    DiscontinuityList& compatible_list)
{
    // Creates a new shock based on the left-most and right-most compatible
    auto new_shock = sh.create_shock(compatible_list[0]->left(),
        compatible_list.back()->right(), compatible_list[0]->type);

    // New position is the mean of all compatible
    double mean_frac = 0.0;
    for (auto& disc : compatible_list) {
        mean_frac += disc->frac;
    }
    mean_frac /= compatible_list.size();
    new_shock.frac = mean_frac;

    // Index and theta are the same as all
    new_shock.ind = compatible_list[0]->ind;
    new_shock.theta = compatible_list[0]->theta;

    // Disconnect entries in shock_point_c
    for (auto& it : compatible_list) {
        auto* sp = dynamic_cast<ShockDiscontinuity*>(it);
        sp->is_connected = false;
    }

    sh.prepare_values(new_shock);
    temp_shock_c.push_back(new_shock);
}

void ShockGrid::find_shocks_to_merge(DiscontinuityList* low_p_left,
    DiscontinuityList* low_p_right, DiscontinuityList& disc_list)
{
    low_p_left->clear();
    low_p_right->clear();
    for (auto& disc : disc_list) {
        auto* sp = dynamic_cast<ShockDiscontinuity*>(disc);
        if (sp != nullptr) {
            if (sp->low_pressure_side == 'l') {
                low_p_left->push_back(sp);
            }
            if (sp->low_pressure_side == 'r') {
                low_p_right->push_back(sp);
            }
        }
    }
}

void ShockGrid::update_values(ShockGrid* grid_to_update_from)
{
    // Grid values update
    CartesianGrid::update_values(grid_to_update_from);
    // Shocks update
    shock_points_c = grid_to_update_from->shock_points_c;
    // Updates maps
    clear_discontinuity_map();
    ShockGrid::fill_discontinuity_map();
}

void ShockGrid::move_shocks(double dt)
{
    for (auto& sp : shock_points_c) {
        double d_frac;
        if (sp.is_x()) {
            d_frac = sp.w * dt / (cos(sp.theta) * dx);
        }
        else {
            d_frac = sp.w * dt / (sin(sp.theta) * dy);
        }
        sp.frac += d_frac;
    }
}

bool ShockGrid::remove_weak_shocks()
{
    const unsigned int orig_vec_size = shock_points_c.size();

    shock_points_c.erase(
        std::remove_if(shock_points_c.begin(), shock_points_c.end(),
            [](ShockDiscontinuity& shock) { return shock.is_weak(); }),
        shock_points_c.end());
    // If a shock is removed we must update the maps
    return (shock_points_c.size() != orig_vec_size);
}

bool ShockGrid::remove_disconnected_shocks()
{
    const unsigned int orig_vec_size = shock_points_c.size();

    shock_points_c.erase(
        std::remove_if(shock_points_c.begin(), shock_points_c.end(),
            [](ShockDiscontinuity& shock) { return !shock.is_connected; }),
        shock_points_c.end());
    // If a shock is removed we must update the maps
    return (shock_points_c.size() != orig_vec_size);
}

bool ShockGrid::detect_new_shocks(double dt)
{
    bool reset_map = false;
    // Detect new shocks
    for (int ind = 0; ind < nPointsTotal; ind++) {
        auto ind_x = indJPlusOne(ind);
        reset_map |= create_x_shock(dt, ind, ind_x);

        auto ind_y = indIPlusOne(ind);
        reset_map |= create_y_shock(dt, ind, ind_y);
    }
    // If a new shock is found, reset the map
    return reset_map;
}

bool ShockGrid::create_x_shock(double dt, int ind, int ind_x)
{
    if (ind_x >= 0) {
        if (ind_x < ind) {
            auto aux = ind;
            ind = ind_x;
            ind_x = aux;
        }
        if (!has_discont_x(ind)
            and flag_functions::point_type(flag(ind)) != SOLID_POINT
            and flag_functions::point_type(flag(ind_x)) != SOLID_POINT) {
            auto new_shock
                = sh.create_shock(values(ind), values(ind_x), 'x', dx / dt);
            new_shock.ind = ind;
            new_shock.frac = 0.5;
            if (new_shock.is_strong()) {
                shock_points_c.push_back(new_shock);
                return true;
            }
        }
    }
    return false;
}

bool ShockGrid::create_y_shock(double dt, int ind, int ind_y)
{
    if (ind_y >= 0) {
        if (ind_y < ind) {
            auto aux = ind;
            ind = ind_y;
            ind_y = aux;
        }
        if (!has_discont_y(ind)
            and flag_functions::point_type(flag(ind)) != SOLID_POINT
            and flag_functions::point_type(flag(ind_y)) != SOLID_POINT) {
            auto new_shock
                = sh.create_shock(values(ind), values(ind_y), 'y', dy / dt);
            new_shock.ind = ind;
            new_shock.frac = 0.5;
            if (new_shock.is_strong()) {
                shock_points_c.push_back(new_shock);
                return true;
            }
        }
    }
    return false;
}

bool ShockGrid::detect_new_connecting_shocks(
    double dt, std::vector<int>& check_x_dir, std::vector<int>& check_y_dir)
{
    bool connected = false;
    /*2*dt makes it detect weaker shocks*/
    for (auto& ind : check_x_dir) {
        connected |= create_x_shock(2 * dt, ind, indJPlusOne(ind));
        connected |= create_x_shock(2 * dt, ind, indJMinusOne(ind));
    }

    for (auto& ind : check_y_dir) {
        connected |= create_y_shock(2 * dt, ind, indIPlusOne(ind));
        connected |= create_y_shock(2 * dt, ind, indIMinusOne(ind));
    }
    return connected;
}

void ShockGrid::compute_shock_angles()
{
    for (auto& sp : shock_points_c) {
        if (sp.type == 'x') {
            compute_shock_theta_x(&sp);
        }
        else {
            compute_shock_theta_y(&sp);
        }
    }
}

void ShockGrid::compute_shock_theta_x(ShockDiscontinuity* sp)
{
    std::vector<ShockDiscontinuity> compatible_up;
    std::vector<ShockDiscontinuity> compatible_down;

    find_compatible_shocks_x(sp, &compatible_up, &compatible_down);

    double mean_pos_x_up, mean_pos_y_up;
    double mean_pos_x_down, mean_pos_y_down;
    if (compatible_up.empty() and compatible_down.empty()) {
        sp->is_connected = false;
        return;
    }
    if (compatible_up.empty()) {
        mean_pos_x_up = disc_X(sp);
        mean_pos_y_up = disc_Y(sp);
        auto mean_down = compute_mean_position(compatible_down);
        mean_pos_x_down = mean_down.first;
        mean_pos_y_down = mean_down.second;
    }
    else if (compatible_down.empty()) {
        mean_pos_x_down = disc_X(sp);
        mean_pos_y_down = disc_Y(sp);
        auto mean_up = compute_mean_position(compatible_up);
        mean_pos_x_up = mean_up.first;
        mean_pos_y_up = mean_up.second;
    }
    else {
        auto mean_down = compute_mean_position(compatible_down);
        mean_pos_x_down = mean_down.first;
        mean_pos_y_down = mean_down.second;
        auto mean_up = compute_mean_position(compatible_up);
        mean_pos_x_up = mean_up.first;
        mean_pos_y_up = mean_up.second;
    }
    auto delta_x = mean_pos_x_up - mean_pos_x_down;
    auto delta_y = mean_pos_y_up - mean_pos_y_down;
    /*atan2(delta_y,delta_x) gives the slope of the line along the shock,
     * but we need the normal, so compute atan2(-delta_x,delta_y) instead*/
    auto normal_theta = atan2(-delta_x, delta_y);
    sp->theta = normal_theta;
    sh.fix_theta(*sp);
}

void ShockGrid::find_compatible_shocks_x(ShockDiscontinuity* sp,
    std::vector<ShockDiscontinuity>* compatible_up,
    std::vector<ShockDiscontinuity>* compatible_down)
{
    auto* discs_x = discontinuity_map_x();
    auto* discs_y = discontinuity_map_y();
    /*
     *   X__b__X
     *   |     |
     *   a     c
     *   |     |
     *   X__O__X
     *   |     |
     *   d     f
     *   |     |
     *   X__e__X
     */

    // Position a
    if (has_discont_y(sp->ind)) {
        extend_compatible_discs(
            compatible_up, sp->low_pressure_side, sp->ind, discs_y);
    }
    // Position b
    auto ind_b = indIPlusOne(sp->ind);
    if (has_discont_x(ind_b)) {
        extend_compatible_discs(
            compatible_up, sp->low_pressure_side, ind_b, discs_x);
    }
    // Position c
    auto ind_c = indJPlusOne(sp->ind);
    if (has_discont_y(ind_c)) {
        auto comp_side = (sp->low_pressure_side == 'l') ? 'r' : 'l';
        extend_compatible_discs(compatible_up, comp_side, ind_c, discs_y);
    }

    auto i_minus = indIMinusOne(sp->ind);
    // Position d
    if (has_discont_y(i_minus)) {
        auto comp_side = (sp->low_pressure_side == 'l') ? 'r' : 'l';
        extend_compatible_discs(compatible_down, comp_side, i_minus, discs_y);
    }
    // Position e
    if (has_discont_x(i_minus)) {
        extend_compatible_discs(
            compatible_down, sp->low_pressure_side, i_minus, discs_x);
    }
    // Position f
    auto ind_f = indJPlusOne(i_minus);
    if (has_discont_y(ind_f)) {
        extend_compatible_discs(
            compatible_down, sp->low_pressure_side, ind_f, discs_y);
    }
}

void ShockGrid::compute_shock_theta_y(ShockDiscontinuity* sp)
{
    std::vector<ShockDiscontinuity> compatible_left;
    std::vector<ShockDiscontinuity> compatible_right;

    find_compatible_shocks_y(sp, &compatible_left, &compatible_right);

    double mean_pos_x_left, mean_pos_y_left;
    double mean_pos_x_right, mean_pos_y_right;
    if (compatible_left.empty() and compatible_right.empty()) {
        sp->is_connected = false;
        return;
    }
    if (compatible_left.empty()) {
        mean_pos_x_left = disc_X(sp);
        mean_pos_y_left = disc_Y(sp);
        auto mean_right = compute_mean_position(compatible_right);
        mean_pos_x_right = mean_right.first;
        mean_pos_y_right = mean_right.second;
    }
    else if (compatible_right.empty()) {
        mean_pos_x_right = disc_X(sp);
        mean_pos_y_right = disc_Y(sp);
        auto mean_left = compute_mean_position(compatible_left);
        mean_pos_x_left = mean_left.first;
        mean_pos_y_left = mean_left.second;
    }
    else {
        auto mean_right = compute_mean_position(compatible_right);
        mean_pos_x_right = mean_right.first;
        mean_pos_y_right = mean_right.second;
        auto mean_left = compute_mean_position(compatible_left);
        mean_pos_x_left = mean_left.first;
        mean_pos_y_left = mean_left.second;
    }
    auto delta_x = mean_pos_x_left - mean_pos_x_right;
    auto delta_y = mean_pos_y_left - mean_pos_y_right;
    /*atan2(delta_y,delta_x) gives the slope of the line along the shock,
     * but we need the normal, so compute atan2(-delta_x,delta_y) instead*/
    auto normal_theta = atan2(-delta_x, delta_y);
    sp->theta = normal_theta;
    sh.fix_theta(*sp);
}

void ShockGrid::find_compatible_shocks_y(ShockDiscontinuity* sp,
    std::vector<ShockDiscontinuity>* compatible_left,
    std::vector<ShockDiscontinuity>* compatible_right)
{
    auto* discs_x = discontinuity_map_x();
    auto* discs_y = discontinuity_map_y();
    /*
     *   X__c__X__f__X
     *   |     |     |
     *   b     O     e
     *   |     |     |
     *   X__a__X__d__X
     */

    // Position a
    auto ind_a = indJMinusOne(sp->ind);
    if (has_discont_x(ind_a)) {
        auto comp_side = (sp->low_pressure_side == 'l') ? 'r' : 'l';
        extend_compatible_discs(compatible_left, comp_side, ind_a, discs_x);
    }
    // Position b
    auto ind_b = indJMinusOne(sp->ind);
    if (has_discont_y(ind_b)) {
        extend_compatible_discs(
            compatible_left, sp->low_pressure_side, ind_b, discs_y);
    }
    // Position c
    auto ind_c = indIPlusOne(ind_b);
    if (has_discont_x(ind_c)) {
        extend_compatible_discs(
            compatible_left, sp->low_pressure_side, ind_c, discs_x);
    }

    // Position d
    auto ind_d = sp->ind;
    if (has_discont_x(ind_d)) {
        extend_compatible_discs(
            compatible_right, sp->low_pressure_side, ind_d, discs_x);
    }
    // Position e
    auto ind_e = indJPlusOne(sp->ind);
    if (has_discont_y(ind_e)) {
        extend_compatible_discs(
            compatible_right, sp->low_pressure_side, ind_e, discs_y);
    }
    // Position f
    auto ind_f = indIPlusOne(sp->ind);
    if (has_discont_x(ind_f)) {
        auto comp_side = (sp->low_pressure_side == 'l') ? 'r' : 'l';
        extend_compatible_discs(compatible_right, comp_side, ind_f, discs_x);
    }
}

std::pair<double, double> ShockGrid::compute_mean_position(
    const std::vector<ShockDiscontinuity>& comp_discs)
{
    double mean_x = 0, mean_y = 0;
    for (auto& comp : comp_discs) {
        mean_x += disc_X(&comp);
        mean_y += disc_Y(&comp);
    }
    return std::make_pair(
        mean_x / comp_discs.size(), mean_y / comp_discs.size());
}

void ShockGrid::extend_compatible_discs(
    std::vector<ShockDiscontinuity>* compatible_discs, char compatible_side,
    int ind, DiscontinuityMap* disc_map)
{
    for (auto* dp : (*disc_map)[ind]) {
        auto* sp = dynamic_cast<ShockDiscontinuity*>(dp);
        if (sp != nullptr) {
            if (sp->low_pressure_side == compatible_side) {
                compatible_discs->push_back(*sp);
            }
        }
    }
}

void ShockGrid::extrapolate_to_shocks()
{
    // Updates all shock points with extrapolated values from the fluid from
    // both sides
    for (auto& sp : shock_points_c) {
        DiscontinuityMap* dir_map;
        int shift_minus;
        int shift_plus;
        int shift_plus_plus;

        std::tie(dir_map, shift_minus, shift_plus, shift_plus_plus)
            = get_map_and_shifts(sp);
        DiscontinuityList& local_disc_list = (*dir_map)[sp.ind];
        if (&sp == local_disc_list[0]) {
            extrapolate_points_left(sp, shift_minus, dir_map);
        }
        if (&sp == local_disc_list.back()) {
            extrapolate_points_right(sp, shift_plus, shift_plus_plus, dir_map);
        }
    }
}

void ShockGrid::extrapolate_points_left(
    ShockDiscontinuity& sp, int shift_minus, DiscontinuityMap* dir_map)
{
    if (dir_map->find(shift_minus) == dir_map->end()) {
        sp.left() = values(sp.ind)
            + (values(sp.ind) - values(shift_minus)) * sp.frac * (1 - sp.frac);
    }
    else {
        sp.left() = values(sp.ind);
    }
}

void ShockGrid::extrapolate_points_right(ShockDiscontinuity& sp, int shift_plus,
    int shift_plus_plus, DiscontinuityMap* dir_map)
{
    if (dir_map->find(shift_plus) == dir_map->end()
        and ind_is_valid(shift_plus)) {
        sp.right() = values(shift_plus)
            + (values(shift_plus) - values(shift_plus_plus)) * (1 - sp.frac)
                * sp.frac;
    }
    else if (ind_is_valid(shift_plus)) {
        sp.right() = values(shift_plus);
    }
    else {
        std::cerr << "Should not be here. An entry was found in the "
                     "discontinuity map at an invalid point!"
                  << std::endl;
    }
}

void ShockGrid::apply_rankine_hugoniot_jump_conditions()
{
    for (auto& sp : shock_points_c) {
        sh.prepare_values(sp);
    }
}

bool ShockGrid::delete_shock_near_wall()
{
    bool shock_deleted = false;
    auto* disc_map = discontinuity_map_x();
    for (auto& disc_x : (*disc_map)) {
        shock_deleted |= find_wall_shock(disc_x.second);
    }
    disc_map = discontinuity_map_y();
    for (auto& disc_y : (*disc_map)) {
        shock_deleted |= find_wall_shock(disc_y.second);
    }
    return shock_deleted;
}

bool ShockGrid::find_wall_shock(DiscontinuityList& disc_list)
{
    bool has_wall = false;
    bool shock_deleted = false;

    for (auto& disc : disc_list) {
        auto* wall = dynamic_cast<BodyDiscontinuity*>(disc);
        if (wall != nullptr) {
            has_wall = true;
        }
    }
    if (has_wall) {
        for (auto& disc : disc_list) {
            auto* shock = dynamic_cast<ShockDiscontinuity*>(disc);
            if (shock != nullptr) {
                shock->is_connected = false;
                shock_deleted = true;
            }
        }
    }

    return shock_deleted;
}

void ShockGrid::specific_print()
{
    std::string number = std::to_string(counter);
    std::string padded_number = std::string(8 - number.length(), '0') + number;
    auto file_name = base_path + "shock" + "_" + padded_number + ".txt";
    std::ofstream output(file_name);
    output << "x,y" << std::endl;
    for (auto& disc : shock_points_c) {
        output << disc_X(&disc) << "," << disc_Y(&disc) << std::endl;
    }
    counter++;
}
