/*!
 * \file derivatives_factory.hpp
 *
 * \brief Headers for create_derivative()
 *
 */
#ifndef DERIVATIVES_FACTORY_HPP
#define DERIVATIVES_FACTORY_HPP

class Derivatives;
class CartesianGrid;
struct PointFunctions;

#include <memory>
#include <string>

/**
 * @brief Creates a pointer to concrete implementation of a derivative
 *
 * @param der_type Regular or Karagiozis
 * @param pf Used to de-reference pointer-to-member-function
 * @param grid Used to get shiftX(Y) and dx(y)
 * @param order Maximum order of accuracy
 *
 * @return shared_ptr to a derivative object implementation
 */
std::shared_ptr<Derivatives> create_derivative(const std::string& der_type,
    PointFunctions& pf, CartesianGrid& grid, const int order = 2);

#endif /* DERIVATIVES_FACTORY_HPP */
