#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXInterpolant.h"

#include <cmath>

#include "CubicInterpolation/Axis.h"

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

double CrossSectionDNDXInterpolant::Calculate(double energy)
{
    return Calculate( energy, 1.0 );
}

double CrossSectionDNDXInterpolant::Calculate(double energy, double v)
{
    auto lim = GetIntegrationLimits(energy);
    v = retransform_relativ_loss(lim[MIN], lim[MAX], v);
    return interpolant.evaluate(std::array<double, 2> { energy, v });
}

double CrossSectionDNDXInterpolant::GetUpperLimit(double energy, double rate)
{
    auto lim = GetIntegrationLimits(energy);
    auto guess = std::array<double, 2> { energy, 0. };
    auto v = cubic_splines::find_parameter(interpolant, rate, guess, 1);
    return transform_relativ_loss(lim[MIN], lim[MAX], v);
}
