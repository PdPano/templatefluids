#ifndef MINIMAL_FILTER_HPP
#define MINIMAL_FILTER_HPP

#include <vector>

class CartesianGrid;
class MinimalFilter {
public:
    MinimalFilter() = default;
    virtual void filter_line(const std::vector<double>& input,
        std::vector<double>* filtered, const std::vector<int>& flags,
        const int len, const int step = 1)
        = 0;
    virtual void filter_grid(CartesianGrid* grid) = 0;
};

#endif /* MINIMAL_FILTER_HPP */
