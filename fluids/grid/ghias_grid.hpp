/**
 * \file ghias_grid.hpp
 * @brief Header for GhiasGrid class
 */
#ifndef GHIAS_GRID_HPP
#define GHIAS_GRID_HPP

#include "../utils/ghias_ghost_point_def.hpp"
#include "../utils/point_functions.hpp"
#include "cartesian_grid.hpp"

// Function from LAPACK
extern "C" void dgesv_(int* lines, int* nrhs, double* mat, int* lda, int* ipiv,
    double* sol, int* ldb, int* info);

/**
 * \class GhiasGrid
 * @brief Implements the ghost point method for handling immersed interfaces
 */
class GhiasGrid : public CartesianGrid {
public:
    GhiasGrid(Options& opt);
    GhiasGrid(Reader reader, Options& opt);
    virtual void grid_specific_update() override;

private:
    /**
     * Container with information relevant to compute the ghost point values
     */
    std::vector<GhiasGhostPoint> ghias_points_c;
    PointFunctions pf;
    /**
     * @name Interpolation functions
     * @param gp Current ghost point
     * @param matrix Matrix of the linear system
     * @param sol Vector of the linear system. Holds the interpolation
     * coefficients after system solution
     * @param line Line to write it
     * @{ */

    /**
     * @brief Generates the linear system matrix
     *
     * @param gp Current ghost point
     * @param matrix Where to write it
     * @param type Condition type ('d' for Dirichlet; 'n' for Neumann)
     */
    void fill_interpolation_matrix(
        GhiasGhostPoint& gp, std::vector<double>& matrix, char type);

    /**
     * @brief Each line in the matrix is filled individualy. This handles one
     * line for the Dirichlet case
     */
    void fill_matrix_line_dirichlet(
        GhiasGhostPoint& gp, std::vector<double>& matrix, int line);
    /**
     * @brief Each line in the matrix is filled individualy. This handles one
     * line for the Neumann case
     */
    void fill_matrix_line_neumann(
        GhiasGhostPoint& gp, std::vector<double>& matrix, int line);

    /**
     * @brief Generates the linear system vector
     *
     * @param type Condition type ('d' for Dirichlet, 'n' for Neumann)
     */
    void fill_interpolation_vector(
        GhiasGhostPoint& gp, std::vector<double>& sol, char type);

    /**
     * @brief Each line is filled individualy. This handles one line for the
     * Neumann case (rho and T have zero normal derivative)
     *
     */
    void fill_interpolation_vector_rho_T(
        GhiasGhostPoint& gp, std::vector<double>& sol, int line);
    /**
     * @brief Each line is filled individualy. This handles one line for the
     * Dirichlet case (ru and rv are zero at the boundary)
     */
    void fill_interpolation_vector_ru_rv(
        GhiasGhostPoint& gp, std::vector<double>& sol, int line);

    /**
     * @brief Calls the LAPACK linear system solver
     */
    void compute_interpolation(
        std::vector<double>& matrix, std::vector<double>& sol);

    /**
     * @brief Compute ghost point values
     *
     * The solution to the system gives the coefficients of the interpolating
     * function. This applies the computed values at the correct point
     *
     * @return pair of doubles. Either (rho,T) or (ru,rv)
     */
    std::pair<double, double> compute_ghost_point_values_from_interpolation(
        GhiasGhostPoint& gp, std::vector<double>& sol);

    /**
     * @brief Transform fluid values in grid values for Dirichlet case
     */
    void compute_dirichlet_values(GhiasGhostPoint& gp,
        std::vector<double>& matrix, std::vector<double>& sol, Point& res);
    /**
     * @brief Transform fluid values in grid values for Neumann case
     */
    void compute_neumann_values(GhiasGhostPoint& gp,
        std::vector<double>& matrix, std::vector<double>& sol, Point& res);
    /**  @} */
};

#endif /* GHIAS_GRID_HPP */
