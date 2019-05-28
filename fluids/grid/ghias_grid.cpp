/**
 * \file ghias_grid.cpp
 * @brief Implementation of GhiasGrid class
 */
#include "ghias_grid.hpp"
#include "../input_output/options.hpp"
#include "../input_output/readers/reader.hpp"

GhiasGrid::GhiasGrid(Options& opt)
    : GhiasGrid(Reader(opt), opt)
{
}

GhiasGrid::GhiasGrid(Reader reader, Options& opt)
    : CartesianGrid(reader)
    , ghias_points_c(reader.ghias_ghost_points())
    , pf(opt.mach(), opt.gam())
{
}

void GhiasGrid::grid_specific_update()
{
    // This assumes that (rho, T) use Neumann conditions and (ru, rv) use
    // Dirichlet conditions
    // Changing it is not too hard, but I'm not doing it
    static std::vector<double> matrix(4 * 4);
    static std::vector<double> sol(2 * 4);
    Point res;
    for (auto& gp : ghias_points_c) {
        compute_dirichlet_values(gp, matrix, sol, res);
        compute_neumann_values(gp, matrix, sol, res);
        this->set_values(res, gp.ind);
    }
}

void GhiasGrid::compute_dirichlet_values(GhiasGhostPoint& gp,
    std::vector<double>& matrix, std::vector<double>& sol, Point& res)
{
    fill_interpolation_matrix(gp, matrix, 'd'); // RU and RV are dirichlet
    fill_interpolation_vector(gp, sol, 'd');
    compute_interpolation(matrix, sol);
    auto dir = compute_ghost_point_values_from_interpolation(gp, sol);
    res.set_ru(-dir.first);
    res.set_rv(-dir.second);
}

void GhiasGrid::compute_neumann_values(GhiasGhostPoint& gp,
    std::vector<double>& matrix, std::vector<double>& sol, Point& res)
{
    fill_interpolation_matrix(gp, matrix, 'n'); // rho and T are neumann
    fill_interpolation_vector(gp, sol, 'n');
    compute_interpolation(matrix, sol);
    auto neu = compute_ghost_point_values_from_interpolation(gp, sol);
    res.set_rho(neu.first);
    res.set_e(pf.e_from_T(res, neu.second));
}

void GhiasGrid::fill_interpolation_matrix(
    GhiasGhostPoint& gp, std::vector<double>& matrix, char type)
{
    // Each ghost point has 4 entries listing where the values used during the
    // interpolation comes from. Each point corresponds to one line in the
    // matrix
    for (int p = 0; p < 4; p++) {
        if (gp.is_fluid(p) or type == 'd') {
            fill_matrix_line_dirichlet(gp, matrix, p);
        }
        else {
            fill_matrix_line_neumann(gp, matrix, p);
        }
    }
}

void GhiasGrid::fill_matrix_line_dirichlet(
    GhiasGhostPoint& gp, std::vector<double>& matrix, int line)
{
    // Fortran index here because of lapack
    matrix[line + 4 * 0] = 1.0;
    matrix[line + 4 * 1] = gp.neighbors_x[line]; // NOLINT
    matrix[line + 4 * 2] = gp.neighbors_y[line]; // NOLINT
    matrix[line + 4 * 3]
        = gp.neighbors_x[line] * gp.neighbors_y[line]; // NOLINT
}

void GhiasGrid::fill_matrix_line_neumann(
    GhiasGhostPoint& gp, std::vector<double>& matrix, int line)
{
    // Fortran index here because of lapack
    matrix[line + 4 * 0] = 0.0;
    matrix[line + 4 * 1] = gp.neighbors_nx[line]; // NOLINT
    matrix[line + 4 * 2] = gp.neighbors_ny[line]; // NOLINT
    matrix[line + 4 * 3]
        = gp.neighbors_nx[line] * gp.neighbors_y[line]  // NOLINT
        + gp.neighbors_ny[line] * gp.neighbors_x[line]; // NOLINT
}

// Only allows for Dirichlet=0 and Neumann=0 for now
void GhiasGrid::fill_interpolation_vector(
    GhiasGhostPoint& gp, std::vector<double>& sol, char type)
{
    for (int p = 0; p < 4; p++) {
        if (gp.is_fluid(p)) {
            if (type == 'd') {
                fill_interpolation_vector_ru_rv(gp, sol, p);
            }
            else {
                fill_interpolation_vector_rho_T(gp, sol, p);
            }
        }
        else {
            sol[p + 4 * 0] = 0.0;
            sol[p + 4 * 1] = 0.0;
        }
    }
}

void GhiasGrid::fill_interpolation_vector_rho_T(
    GhiasGhostPoint& gp, std::vector<double>& sol, int line)
{
    auto p = values(gp.neighbors_inds[line]); // NOLINT
    sol[line + 4 * 0] = pf.rho(p);
    sol[line + 4 * 1] = pf.temperature(p);
}

void GhiasGrid::fill_interpolation_vector_ru_rv(
    GhiasGhostPoint& gp, std::vector<double>& sol, int line)
{
    auto p = values(gp.neighbors_inds[line]); // NOLINT
    sol[line + 4 * 0] = pf.ru(p);
    sol[line + 4 * 1] = pf.rv(p);
}

void GhiasGrid::compute_interpolation(
    std::vector<double>& matrix, std::vector<double>& sol)
{
    static std::vector<int> ipiv(4);
    static int dim = 4;
    static int nrhs = 2;
    static int info;
    // Call to LAPACK
    dgesv_(&dim, &nrhs, &*matrix.begin(), &dim, &*ipiv.begin(), &*sol.begin(),
        &dim, &info);
    if (info != 0) {
        std::cerr << "[GHIAS_GRID] Problems with system solution!" << std::endl;
        std::cerr << "lapack info = " << info << std::endl;
    }
}

std::pair<double, double>
GhiasGrid::compute_ghost_point_values_from_interpolation(
    GhiasGhostPoint& gp, std::vector<double>& sol)
{
    auto x = gp.image_coordinate[0];
    auto y = gp.image_coordinate[1];
    auto res = [&](int i) {
        return sol[0 + 4 * i] + sol[1 + 4 * i] * x + sol[2 + 4 * i] * y
            + sol[3 + 4 * i] * x * y;
    };
    return std::make_pair(res(0), res(1));
}
