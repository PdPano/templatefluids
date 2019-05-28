#include "minimal_filter_factory.hpp"

#include "minimal_filter_3_moments.hpp"

#include <iostream>

std::shared_ptr<MinimalFilter> create_minimal_filter(int order)
{
    if (order == 3) {
        return std::make_shared<MinimalFilter3Moments>();
    }
    std::cerr << "Filter of order " << order
              << " not available. Using order=3 instead!!!" << std::endl;
    return std::make_shared<MinimalFilter3Moments>();
}
