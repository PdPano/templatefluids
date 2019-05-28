#include "weno_convection.hpp"
#include "../grid/cartesian_grid.hpp"
#include "../utils/flag_handler.hpp"
#include "../utils/flux_def.hpp"
#include "../utils/operators_overloads.hpp"
#include "../utils/point_def.hpp"
#include "../utils/point_functions.hpp"
#include "flux_functions/flux_interface.hpp"

#include <numeric>
#include <utility>

WenoConvection::WenoConvection(PointFunctions& pf_in,
    std::shared_ptr<FluxInterface> flux_in,
    std::shared_ptr<Convection> backup_convection_in)
    : Convection(pf_in, nullptr)
    , flux(std::move(flux_in))
    , backup_convection(std::move(backup_convection_in))
{
}

Flux WenoConvection::convection_x(const CartesianGrid& grid, int ind) const
{
    auto der_flag = flag_functions::derx(grid.flag(ind));
    if (std::min(
            flag_functions::left(der_flag), flag_functions::right(der_flag))
        >= 3) {
        auto shiftX = grid.shiftX();
        auto f_plus_i
            = flux_at_i_half(grid, ind, &FluxInterface::fluxXPositive, shiftX);
        auto f_plus_i_m = flux_at_i_half(
            grid, ind - shiftX, &FluxInterface::fluxXPositive, shiftX);
        auto f_minus_i
            = flux_at_i_half(grid, ind, &FluxInterface::fluxXNegative, -shiftX);
        auto f_minus_i_p = flux_at_i_half(
            grid, ind + shiftX, &FluxInterface::fluxXNegative, -shiftX);

        return -((f_plus_i - f_plus_i_m) + (f_minus_i_p - f_minus_i)) / grid.dx;
    }
    return backup_convection->convection_x(grid, ind);
}

Flux WenoConvection::convection_y(const CartesianGrid& grid, int ind) const
{
    auto der_flag = flag_functions::dery(grid.flag(ind));
    if (std::min(
            flag_functions::left(der_flag), flag_functions::right(der_flag))
        >= 3) {
        auto shiftY = grid.shiftY();
        auto f_plus_i
            = flux_at_i_half(grid, ind, &FluxInterface::fluxYPositive, shiftY);
        auto f_plus_i_m = flux_at_i_half(
            grid, ind - shiftY, &FluxInterface::fluxYPositive, shiftY);
        auto f_minus_i
            = flux_at_i_half(grid, ind, &FluxInterface::fluxYNegative, -shiftY);
        auto f_minus_i_p = flux_at_i_half(
            grid, ind + shiftY, &FluxInterface::fluxYNegative, -shiftY);

        return -((f_plus_i - f_plus_i_m) + (f_minus_i_p - f_minus_i)) / grid.dy;
    }
    return backup_convection->convection_y(grid, ind);
}

WenoConvection::BetasArray WenoConvection::compute_betas(
    const ComponentsArray& vals) const
{
    // This is a convenient offset to keep the notation consistent with
    // the WENO papers
    const int i = 2;
    BetasArray betas{};
    for (int comp = 0; comp < 4; comp++) {
        auto& f = vals.at(comp);
        auto& b = betas.at(comp);
        b[0] = 13 / 12. * pow(f[i - 2] - 2 * f[i - 1] + f[i], 2)
            + 1 / 4. * pow(f[i - 2] - 4 * f[i - 1] + 3 * f[i], 2);
        b[1] = 13 / 12. * pow(f[i - 1] - 2 * f[i] + f[i + 1], 2)
            + 1 / 4. * pow(f[i - 1] - f[i + 1], 2);
        b[2] = 13 / 12. * pow(f[i] - 2 * f[i + 1] + f[i + 2], 2)
            + 1 / 4. * pow(f[i + 2] - 4 * f[i + 1] + 3 * f[i], 2);
    }
    return betas;
}

WenoConvection::SigmasArray WenoConvection::compute_sigmas(
    const BetasArray& betas) const
{
    const double eps = 1e-30;
    SigmasArray sigmas = {};
    for (int i = 0; i < 4; i++) {
        auto& s = sigmas.at(i);
        auto& b = betas.at(i);
        s[0] = (1 / 10.) / pow(eps + b[0], 2);
        s[1] = (3 / 5.) / pow(eps + b[1], 2);
        s[2] = (3 / 10.) / pow(eps + b[2], 2);
    }

    return sigmas;
}

WenoConvection::OmegasArray WenoConvection::compute_omegas(
    const SigmasArray& sigmas) const
{
    OmegasArray omegas = {};
    for (int i = 0; i < 4; i++) {
        auto& w = omegas.at(i);
        auto& s = sigmas.at(i);
        const double sum = s[0] + s[1] + s[2];
        w[0] = s[0] / sum;
        w[1] = s[1] / sum;
        w[2] = s[2] / sum;
    }

    return omegas;
}

WenoConvection::ApproximationArray WenoConvection::compute_approximations(
    const ComponentsArray& vals) const
{
    const int i = 2;
    ApproximationArray approx = {};
    for (int comp = 0; comp < 4; comp++) {
        auto& f = vals.at(comp);
        auto& a = approx.at(comp);
        a[0] = 1 / 3. * f[i - 2] - 7 / 6. * f[i - 1] + 11 / 6. * f[i];
        a[1] = -1 / 6. * f[i - 1] + 5 / 6. * f[i] + 1 / 3. * f[i + 1];
        a[2] = 1 / 3. * f[i] + 5 / 6. * f[i + 1] - 1 / 6. * f[i + 2];
    }
    return approx;
}

Flux WenoConvection::flux_at_i_half(const CartesianGrid& grid, int ind,
    WenoConvection::fluxFunc func, int shift) const
{
    Flux res = {0., 0., 0., 0.};
    auto components = build_components_array(grid, ind, func, shift);
    auto sigmas = compute_sigmas(compute_betas(components));
    auto omegas = compute_omegas(sigmas);
    auto approx = compute_approximations(components);

    res.rho = std::inner_product(
        approx[0].begin(), approx[0].end(), omegas[0].begin(), 0.0);
    res.ru = std::inner_product(
        approx[1].begin(), approx[1].end(), omegas[1].begin(), 0.0);
    res.rv = std::inner_product(
        approx[2].begin(), approx[2].end(), omegas[2].begin(), 0.0);
    res.e = std::inner_product(
        approx[3].begin(), approx[3].end(), omegas[3].begin(), 0.0);

    return res;
}

WenoConvection::ComponentsArray WenoConvection::build_components_array(
    const CartesianGrid& grid, int ind, fluxFunc func, int shift) const
{
    auto flux_at = [&](int i) { return ((flux.get())->*func)(grid.values(i)); };
    ComponentsArray components{};
    /* ind - 2*shift, ind-shift,...,ind+2*shift  */
    for (int i = 0; i < 5; i++) {
        auto f = flux_at(ind + (i - 2) * shift);
        components[0].at(i) = f.rho;
        components[1].at(i) = f.ru;
        components[2].at(i) = f.rv;
        components[3].at(i) = f.e;
    }
    return components;
}
