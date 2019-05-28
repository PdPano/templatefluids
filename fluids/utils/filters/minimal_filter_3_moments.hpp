#ifndef MINIMAL_FILTER_3_MOMENTS_HPP
#define MINIMAL_FILTER_3_MOMENTS_HPP

#include "minimal_filter.hpp"

class MinimalFilter3Moments : public MinimalFilter {
public:
    MinimalFilter3Moments() = default;
    void filter_line(const std::vector<double>& input,
        std::vector<double>* filtered, const std::vector<int>& flags,
        const int len, const int step = 1);
    void filter_grid(CartesianGrid* grid);

private:
    using Flags = const std::vector<int>;

    template <typename T>
    T center_filter(int ind, const std::vector<T>& input, const int step);

    template <typename T>
    T left_filter(int ind, const std::vector<T>& input, const int step);

    template <typename T>
    T almost_left_filter(int ind, const std::vector<T>& input, const int step);

    template <typename T>
    T right_filter(
        int ind, const std::vector<T>& input, const int step, int sign = 1);

    template <typename T>
    T almost_right_filter(
        int ind, const std::vector<T>& input, const int step, int sign = 1);

    template <typename T>
    T select_filter(
        int ind, const std::vector<T>& input, Flags& flags, const int step);
};

#endif /* MINIMAL_FILTER_3_MOMENTS_HPP */
