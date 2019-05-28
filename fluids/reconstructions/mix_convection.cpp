#include "mix_convection.hpp"
#include "../derivatives/derivatives.hpp"
#include "../grid/cartesian_grid.hpp"
#include "../utils/flux_def.hpp"
#include "../utils/operators_overloads.hpp"
#include "../utils/point_def.hpp"
#include "../utils/point_functions.hpp"
#include "flux_functions/flux_interface.hpp"

#include <utility>

MixConvection::MixConvection(PointFunctions& pf_in,
    std::shared_ptr<Derivatives> der_in, std::shared_ptr<Convection> main_in,
    std::shared_ptr<Convection> aux_in, double mix_param_in)
    : Convection(pf_in, std::move(der_in))
    , main(std::move(main_in))
    , aux(std::move(aux_in))
    , mix_param(mix_param_in)
{
}

Flux MixConvection::convection_x(const CartesianGrid& grid, int ind) const
{
    return (1 - mix_param) * main->convection_x(grid, ind)
        + mix_param * aux->convection_x(grid, ind);
}

Flux MixConvection::convection_y(const CartesianGrid& grid, int ind) const
{
    return (1 - mix_param) * main->convection_y(grid, ind)
        + mix_param * aux->convection_y(grid, ind);
}
