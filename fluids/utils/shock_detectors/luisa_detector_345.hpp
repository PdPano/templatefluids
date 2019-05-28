#ifndef LUISA_DETECTOR_345_HPP
#define LUISA_DETECTOR_345_HPP

#include "luisa_shock_detector.hpp"
#include "../filters/minimal_filter.hpp"
#include <vector>
#include <memory>

/**
 * \class LuisaDetector345
 * @brief Implements Luisa's shock detector using third, fourth and fifth
 * derivatives
 */

class LuisaDetector345 : public LuisaDetector {
public:
    LuisaDetector345(double sensitivity_in, int nPointsI_in, int nPointsJ_in);
    void detect_shocks(const CartesianGrid& grid);
    const std::set<int>& set_of_shocked_points(void) const
    {
        return shocked_points;
    }

private:
    const double sensitivity;
    const int nPointsI;
    const int nPointsJ;
    std::vector<double> rho_to_derive;
    std::vector<double> rho_to_derive_two;
    std::vector<double> rho_to_load;
    std::vector<double> numerator;
    std::vector<double> denominator;
    std::vector<int> der_flags;
    std::set<int> shocked_points;
    std::shared_ptr<MinimalFilter> minimal_filter;
    void load_data(int line_or_col, char dir, const CartesianGrid& grid);
    void find_shocks(int line_or_col, char dir, std::set<int>* shocked_points,
        const CartesianGrid& grid);
    void derive_vector(std::vector<double>* in, int len, int step);
    void compute_numerator(int len);
    void compute_denominator(int len);
    double generic_first_der(int ind, const std::vector<double>& rho, int step);
};

#endif /* LUISA_DETECTOR_345_HPP */
