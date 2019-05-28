/*!
 * \file derivatives_factory.cpp
 *
 * \brief Implementation of create_derivative()
 *
 * Generates a derivative object of required type
 *
 */
#include "derivatives_factory.hpp"
#include "../grid/cartesian_grid.hpp"
#include "../utils/point_functions.hpp"
#include "../utils/useful_alias.hpp"
#include "irregular_derivatives.hpp"
#include "regular_derivatives.hpp"
#include "regular_derivatives_fourth.hpp"
#include <iostream>
std::shared_ptr<Derivatives> create_derivative(const std::string& der_type,
    PointFunctions& pf, CartesianGrid& grid, const int order)
{

    if (der_type == "REGULAR") {
        if (order == 4) {
            return std::make_shared<RegularDerivativesFourth>(
                pf, grid.shiftX(), grid.shiftY(), grid.dx, grid.dy);
        }
        return std::make_shared<RegularDerivatives>(
            pf, grid.shiftX(), grid.shiftY(), grid.dx, grid.dy);
    }
    if (der_type == "KARAGIOZIS") {
        return std::make_shared<IrregularDerivatives>(
            pf, grid.shiftX(), grid.shiftY(), grid.dx, grid.dy);
    }
    std::cerr << "[DEV-ERROR] Could not create derivatives of type " << der_type
              << std::endl;
    return nullptr;
}
