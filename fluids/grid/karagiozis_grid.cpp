/**
 * \file karagiozis_grid.cpp
 * @brief Implementation of KaragiozisGrid class
 */
#include "karagiozis_grid.hpp"
#include "../input_output/options.hpp"
#include "../input_output/readers/reader.hpp"
#include <algorithm>
#include <iostream>
#include <tuple>
#include <utility>

KaragiozisGrid::KaragiozisGrid(Options& opt)
    : KaragiozisGrid(Reader(opt), opt)
{
}

KaragiozisGrid::KaragiozisGrid(Reader reader, Options& opt)
    : CartesianGrid(reader)
    , pf(opt.mach(), opt.gam())
    , karagiozis_points_c(reader.karagiozis_body_points())
{
    KaragiozisGrid::fill_discontinuity_map();
    KaragiozisGrid::sort_lists();
    std::ofstream outfile("./output/disc_points.txt");
    for (auto& dp : karagiozis_points_c) {
        outfile << disc_X(&dp) << "," << disc_Y(&dp) << std::endl;
    }
}

auto KaragiozisGrid::get_map_and_shifts(BodyDiscontinuity& bd)
{
    int shift_minus;
    int shift_plus;
    int shift_plus_plus;
    if (bd.is_x()) {
        shift_minus = indJMinusOne(bd.ind);
        shift_plus = indJPlusOne(bd.ind);
        shift_plus_plus = indJPlusOne(shift_plus);
        if (!ind_is_valid(shift_minus)) {
            shift_minus = bd.ind;
        }
        if (!ind_is_valid(shift_plus_plus)) {
            shift_plus_plus = shift_plus;
        }
        return std::make_tuple(
            discontinuity_map_x(), shift_minus, shift_plus, shift_plus_plus);
    }
    // Assuming bd is y
    shift_minus = indIMinusOne(bd.ind);
    shift_plus = indIPlusOne(bd.ind);
    shift_plus_plus = indIPlusOne(shift_plus);
    if (!ind_is_valid(shift_minus)) {
        shift_minus = bd.ind;
    }
    if (!ind_is_valid(shift_plus_plus)) {
        shift_plus_plus = shift_plus;
    }
    return std::make_tuple(
        discontinuity_map_y(), shift_minus, shift_plus, shift_plus_plus);
}

void KaragiozisGrid::clear_discontinuity_map()
{
    disc_map_x.clear();
    disc_map_y.clear();
}

void KaragiozisGrid::fill_discontinuity_map()
{
    for (auto& bd : karagiozis_points_c) {
        DiscontinuityMap* dir_map;
        if (bd.is_x()) {
            dir_map = discontinuity_map_x();
        }
        else {
            dir_map = discontinuity_map_y();
        }
        (*dir_map)[bd.ind].push_back(&bd);
        set_points_to_revisit(&bd);
    }
}

void KaragiozisGrid::set_points_to_revisit(GenericDiscontinuity* disc)
{
    int ind = disc->ind;
    int ind_right;
    if (disc->is_x()) {
        /*    O----O
         *    |    |
         *ind O-x--O
         *    |    |
         *    O----O
         */
        ind_right = indJPlusOne(ind);

        safe_insert_to_points_to_revisit(ind);
        safe_insert_to_points_to_revisit(ind_right);
        safe_insert_to_points_to_revisit(indIPlusOne(ind));
        safe_insert_to_points_to_revisit(indIPlusOne(ind_right));
        safe_insert_to_points_to_revisit(indIMinusOne(ind));
        safe_insert_to_points_to_revisit(indIMinusOne(ind_right));
    }
    else if (disc->is_y()) {
        /*   O----O----O
         *   |    x    |
         *   |    |    |
         *   O----O----O
         *       ind
         */
        ind_right = indIPlusOne(ind);

        safe_insert_to_points_to_revisit(ind);
        safe_insert_to_points_to_revisit(ind_right);
        safe_insert_to_points_to_revisit(indJPlusOne(ind));
        safe_insert_to_points_to_revisit(indJPlusOne(ind_right));
        safe_insert_to_points_to_revisit(indJMinusOne(ind));
        safe_insert_to_points_to_revisit(indJMinusOne(ind_right));
    }
}

void KaragiozisGrid::safe_insert_to_points_to_revisit(int ind)
{
    if (ind >= 0) {
        points_to_revisit.insert(ind);
    }
}

void KaragiozisGrid::sort_lists()
{
    auto less_operator = [](GenericDiscontinuity* a, GenericDiscontinuity* b) {
        return (a->frac) < (b->frac);
    };
    for (auto& list : *(discontinuity_map_x())) {
        std::sort(list.second.begin(), list.second.end(), less_operator);
    }
    for (auto& list : *(discontinuity_map_y())) {
        std::sort(list.second.begin(), list.second.end(), less_operator);
    }
}

void KaragiozisGrid::grid_specific_update()
{
    // Updates all body points with extrapolated values from the fluid from both
    // sides
    for (auto& bd : karagiozis_points_c) {
        DiscontinuityMap* dir_map;
        int shift_minus;
        int shift_plus;
        int shift_plus_plus;

        std::tie(dir_map, shift_minus, shift_plus, shift_plus_plus)
            = get_map_and_shifts(bd);
        DiscontinuityList& local_disc_list = (*dir_map)[bd.ind];
        if (&bd == local_disc_list[0]) {
            extrapolate_left(bd, shift_minus, dir_map);
        }
        if (&bd == local_disc_list.back()) {
            extrapolate_right(bd, shift_plus, shift_plus_plus, dir_map);
        }
    }
}

void KaragiozisGrid::extrapolate_left(
    BodyDiscontinuity& bd, int shifted_ind, DiscontinuityMap* dir_map)
{
    double new_rho;
    double new_T = pf.temperature(this->values(bd.ind));
    if (dir_map->find(shifted_ind) == dir_map->end()) {
        new_rho = linear_extrapolation(
            this->rho(shifted_ind), this->rho(bd.ind), bd.frac);
    }
    else {
        new_rho = this->rho(bd.ind);
    }
    bd.left().set_rho(new_rho);
    double new_e = pf.e_from_T(bd.left(), new_T); // Fix the temperature
    bd.left().set_e(new_e);
}

void KaragiozisGrid::extrapolate_right(BodyDiscontinuity& bd, int shift_plus,
    int shift_plus_plus, DiscontinuityMap* dir_map)
{
    double new_rho;
    double new_T = pf.temperature(this->values(shift_plus));
    if (dir_map->find(shift_plus) == dir_map->end()
        and ind_is_valid(shift_plus)) {
        new_rho = linear_extrapolation(
            this->rho(shift_plus_plus), this->rho(shift_plus), 1 - bd.frac);
    }
    else if (ind_is_valid(shift_plus)) {
        new_rho = this->rho(shift_plus);
    }
    else {
        std::cerr << "Should not be here. An entry was found in the "
                     "discontinuity map at an invalid point!"
                  << std::endl;
        new_rho = 1.0;
    }
    bd.right().set_rho(new_rho);
    double new_e = pf.e_from_T(bd.right(), new_T); // Fix the temperature
    bd.right().set_e(new_e);
}

bool KaragiozisGrid::has_discont_x(int ind, int shift) const
{
    auto xmap = discontinuity_map_x();
    auto inds = indIJ(ind);
    auto new_ind = IND(inds.first, inds.second + shift);
    return xmap->find(new_ind) != xmap->end();
}

bool KaragiozisGrid::has_discont_y(int ind, int shift) const
{
    auto ymap = discontinuity_map_y();
    auto inds = indIJ(ind);
    auto new_ind = IND(inds.first + shift, inds.second);
    return ymap->find(new_ind) != ymap->end();
}

// This function is only called after it is known a discontinuity exists at ind
GenericDiscontinuity* KaragiozisGrid::first_disc_x(int ind) const
{
    auto xmap = discontinuity_map_x();
    const DiscontinuityList& vec = xmap->at(ind);
    return vec[0];
}

// This function is only called after it is known a discontinuity exists at ind
GenericDiscontinuity* KaragiozisGrid::last_disc_x(int ind) const
{
    auto xmap = discontinuity_map_x();
    const DiscontinuityList& vec = xmap->at(ind);
    return vec.back();
}

// This function is only called after it is known a discontinuity exists at ind
GenericDiscontinuity* KaragiozisGrid::first_disc_y(int ind) const
{
    auto ymap = discontinuity_map_y();
    const DiscontinuityList& vec = ymap->at(ind);
    return vec[0];
}

// This function is only called after it is known a discontinuity exists at ind
GenericDiscontinuity* KaragiozisGrid::last_disc_y(int ind) const
{
    auto ymap = discontinuity_map_y();
    const DiscontinuityList& vec = ymap->at(ind);
    return vec.back();
}
