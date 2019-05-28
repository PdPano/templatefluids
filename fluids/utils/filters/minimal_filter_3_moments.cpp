#include "minimal_filter_3_moments.hpp"
#include "../../grid/cartesian_grid.hpp"
#include "../flag_handler.hpp"
#include "../operators_overloads.hpp"
#include "../point_def.hpp"

void MinimalFilter3Moments::filter_line(const std::vector<double>& input,
    std::vector<double>* filtered, const std::vector<int>& flags, const int len,
    const int step)
{
    for (int i = 0; i < len; i++) {
        (*filtered)[i] = select_filter<double>(i, input, flags, step);
    }
}

void MinimalFilter3Moments::filter_grid(CartesianGrid* grid)
{
    std::vector<Point> line(std::max(grid->nPointsI, grid->nPointsJ));
    std::vector<int> flags(std::max(grid->nPointsI, grid->nPointsJ));

    /*filter in x direction*/
    for (int i = 0; i < grid->nPointsI; i++) {
        for (int j = 0; j < grid->nPointsJ; j++) {
            int ind = grid->IND(i, j);
            line[j] = grid->values(ind);
            flags[j] = flag_functions::derx(grid->flag(ind));
        }
        for (int j = 0; j < grid->nPointsJ; j++) {
            int ind = grid->IND(i, j);
            grid->set_values(select_filter(j, line, flags, 1), ind);
        }
    }

    /*filter in y direction*/
    for (int j = 0; j < grid->nPointsJ; j++) {
        for (int i = 0; i < grid->nPointsI; i++) {
            int ind = grid->IND(i, j);
            line[i] = grid->values(ind);
            flags[i] = flag_functions::dery(grid->flag(ind));
        }
        for (int i = 0; i < grid->nPointsI; i++) {
            int ind = grid->IND(i, j);
            grid->set_values(select_filter(i, line, flags, 1), ind);
        }
    }
}

template <typename T>
T MinimalFilter3Moments::center_filter(
    int ind, const std::vector<T>& input, const int step)
{
    return -1 / 16. * input[ind - 2 * step] + 1 / 4. * input[ind - 1 * step]
        + 5 / 8. * input[ind] + 1 / 4. * input[ind + 1 * step]
        - 1 / 16. * input[ind + 2 * step];
}

template <typename T>
T MinimalFilter3Moments::left_filter(
    int ind, const std::vector<T>& input, const int step)
{
    return right_filter(ind, input, step, -1);
}

template <typename T>
T MinimalFilter3Moments::almost_left_filter(
    int ind, const std::vector<T>& input, const int step)
{
    return almost_right_filter(ind, input, step, -1);
}

template <typename T>
T MinimalFilter3Moments::right_filter(
    int ind, const std::vector<T>& input, const int step, int sign)
{
    return 15 / 16. * input[ind] + 1 / 4. * input[ind + 1 * sign * step]
        - 3 / 8. * input[ind + 2 * sign * step]
        + 1 / 4. * input[ind + 3 * sign * step]
        - 1 / 16. * input[ind + 4 * sign * step];
}

template <typename T>
T MinimalFilter3Moments::almost_right_filter(
    int ind, const std::vector<T>& input, const int step, int sign)
{
    return 1 / 16. * input[ind - 1 * sign * step] + 3 / 4. * input[ind]
        + 3 / 8. * input[ind + 1 * sign * step]
        - 1 / 4. * input[ind + 2 * sign * step]
        + 1 / 16. * input[ind + 3 * sign * step];
}

template <typename T>
T MinimalFilter3Moments::select_filter(
    int ind, const std::vector<T>& input, Flags& flags, const int step)
{
    auto der_flag = flags[ind];
    if (std::min(
            flag_functions::left(der_flag), flag_functions::right(der_flag))
        > 1) {
        return center_filter(ind, input, step);
    }
    if (flag_functions::left(der_flag) < 1 * step
        and flag_functions::right(der_flag) >= 4 * step) {
        return right_filter(ind, input, step);
    }
    if (flag_functions::left(der_flag) < 2 * step
        and flag_functions::right(der_flag) >= 3 * step) {
        return almost_right_filter(ind, input, step);
    }
    if (flag_functions::right(der_flag) < 1 * step
        and flag_functions::left(der_flag) >= 4 * step) {
        return left_filter(ind, input, step);
    }
    if (flag_functions::right(der_flag) < 2 * step
        and flag_functions::left(der_flag) >= 3 * step) {
        return almost_left_filter(ind, input, step);
    }
    return input[ind];
}
