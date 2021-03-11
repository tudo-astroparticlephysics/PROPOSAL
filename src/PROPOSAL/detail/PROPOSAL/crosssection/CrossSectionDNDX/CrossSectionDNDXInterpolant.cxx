#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXInterpolant.h"

#include <cmath>

#include "CubicInterpolation/Axis.h"
#include "CubicInterpolation/FindParameter.hpp"
using namespace PROPOSAL;

template <>
std::function<double(double, double, double)> transform_loss<ComptonKleinNishina>::func
    = [](double v_cut, double v_max, double v) {
          return transform_loss_linear(v_cut, v_max, v);
      };

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

std::string CrossSectionDNDXInterpolant::gen_path() const
{
    return std::string(InterpolationSettings::TABLES_PATH);
}

std::string CrossSectionDNDXInterpolant::gen_name() const
{
    return std::string("dndx_") + std::to_string(GetHash())
        + std::string(".txt");
}

double CrossSectionDNDXInterpolant::Calculate(double energy)
{
    return Calculate(energy, 1.0);
}

double CrossSectionDNDXInterpolant::Calculate(double energy, double v)
{
    if (energy < lower_energy_lim)
        return 0.;
    auto lim = GetIntegrationLimits(energy);
    v = retransform_relativ_loss(lim.min, lim.max, v);
    return interpolant.evaluate(std::array<double, 2> { energy, v });
}

double CrossSectionDNDXInterpolant::GetUpperLimit(double energy, double rate)
{
    if (energy < lower_energy_lim)
        throw std::invalid_argument("no dNdx for this energy defined.");
    auto lim = GetIntegrationLimits(energy);

    auto initial_guess = cubic_splines::ParameterGuess<std::array<double, 2>> {
        .x = { energy, NAN },
        .n = 1,
    };

    auto v = cubic_splines::find_parameter(interpolant, rate, initial_guess);
    return transform_relativ_loss(lim.min, lim.max, v);
}
