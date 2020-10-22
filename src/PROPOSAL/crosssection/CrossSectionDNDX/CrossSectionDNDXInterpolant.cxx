#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXInterpolant.h"
#include <cmath>

using namespace PROPOSAL;

namespace PROPOSAL {
double transform_relativ_loss(double v_cut, double v_max, double v)
{
    if (v < 0 or v_max == 0)
        return v_cut;
    if (v >= 1)
        return v_max;
    return v_cut * std::exp(v * std::log(v_max / v_cut));
}

double retransform_relativ_loss(double v_cut, double v_max, double v)
{
    if (v <= v_cut)
        return 0;
    if (v >= v_max)
        return 1;
    return std::log(v / v_cut) / std::log(v_max / v_cut);
}
}


unique_ptr<Interpolant> CrossSectionDNDXInterpolant::build_dndx(crosssection::Parametrization const& param, ParticleDef const& p_def)
{
    dNdx_def.xmin = param.GetLowerEnergyLim(p_def);
    dNdx_def.isLog = true;
    dNdx_def.rationalY = true;
    dNdx_def.x2min = 0.0;
    dNdx_def.x2max = 1.0;
    dNdx_def.function2d = [this](double energy, double v) {
        auto integral_lim = GetIntegrationLimits(energy);
        v = transform_relativ_loss(get<MIN>(integral_lim), get<MAX>(integral_lim), v);
        return CrossSectionDNDXIntegral::Calculate(energy, v);
    };

    auto hash_digest = (size_t)0;
    hash_combine(hash_digest, hash_cross_section, dNdx_def.GetHash());

    return Helper::InitializeInterpolation("dNdx", Interpolant2DBuilder(dNdx_def), hash_digest);
}

unique_ptr<Interpolant> CrossSectionDNDXInterpolant::build_dndx1d(crosssection::Parametrization const& param, ParticleDef const& p_def)
{
    dNdx_def.xmin = param.GetLowerEnergyLim(p_def);
    dNdx_def.isLog = true;
    dNdx_def.rationalY = true;
    //TODO: this may be unnecessary as soon as we have a more efficient
    //      implementation of the (2d) interpolation (jm)
    dNdx_def.function1d = [this](double energy) {
        // auto integral_lim = GetIntegrationLimits(energy);
        return dndx->Interpolate(energy, 1.);
    };

    auto hash_digest = (size_t)0;
    hash_combine(hash_digest, hash_cross_section, dNdx_def.GetHash());

    return Helper::InitializeInterpolation("dNdx1d", Interpolant1DBuilder(dNdx_def), hash_cross_section);
}

double CrossSectionDNDXInterpolant::Calculate(double energy)
{
    return dndx1d->Interpolate(energy);
}

double CrossSectionDNDXInterpolant::Calculate(double energy, double v)
{
    auto integral_lim = GetIntegrationLimits(energy);
    v = retransform_relativ_loss(get<MIN>(integral_lim), get<MAX>(integral_lim), v);
    return dndx->Interpolate(energy, v);
}

double CrossSectionDNDXInterpolant::GetUpperLimit(double energy, double rate)
{
    auto integral_lim = GetIntegrationLimits(energy);
    auto v = dndx->FindLimit(energy, rate);
    return transform_relativ_loss(get<MIN>(integral_lim), get<MAX>(integral_lim), v);
}
