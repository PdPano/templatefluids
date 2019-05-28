#ifndef LUISA_SHOCK_DETECTOR_HPP
#define LUISA_SHOCK_DETECTOR_HPP
#include <set>

class CartesianGrid;
class LuisaDetector {
public:
    LuisaDetector() {}
    virtual void detect_shocks(const CartesianGrid& grid) = 0;
    virtual const std::set<int>& set_of_shocked_points() const = 0;
};

#endif /* LUISA_SHOCK_DETECTOR_HPP */
