#include "../input_output/options.hpp"
#include "../utils/point_functions.hpp"
#include "abstract_convection.hpp"
#include "convection_factory.hpp"
#include "flux_functions/flux_factory.hpp"
#include "mix_convection.hpp"
#include "simple_convection.hpp"
#include "simple_flux_convection.hpp"
#include "skew_symmetric.hpp"
#include "split_convection.hpp"
#include "split_convection_cached.hpp"
#include "weno_convection.hpp"

std::shared_ptr<Convection> create_convection(Options& opt, PointFunctions& pf,
    std::shared_ptr<Derivatives> der, std::string overwrite_conv)
{
    if (overwrite_conv == "NONE") {
        overwrite_conv = opt.convection();
    }
    if (overwrite_conv == "SIMPLE") {
        return std::make_shared<SimpleConvection>(pf, der);
    }
    if (overwrite_conv == "SKEW_SYMMETRIC") {
        return std::make_shared<SkewSymmetric>(pf, der);
    }
    if (overwrite_conv == "SPLIT_CONVECTION") {
        return std::make_shared<SplitConvection>(pf, der, create_flux(opt, pf));
    }
    if (overwrite_conv == "SIMPLE_FLUX") {
        return std::make_shared<SimpleFluxConvection>(
            pf, der, create_flux(opt, pf));
    }
    if (overwrite_conv == "WENO_CONVECTION") {
        auto bkp_conv
            = std::make_shared<SplitConvection>(pf, der, create_flux(opt, pf));
        return std::make_shared<WenoConvection>(
            pf, create_flux(opt, pf), bkp_conv);
    }
    if (overwrite_conv == "MIX_CONVECTION") {
        auto main_conv
            = create_convection(opt, pf, der, opt.mix_convection_main());
        auto aux_conv
            = create_convection(opt, pf, der, opt.mix_convection_aux());
        return std::make_shared<MixConvection>(
            pf, der, main_conv, aux_conv, opt.mix_param());
    }
    return nullptr;
}
