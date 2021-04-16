#define CROSSSECTIONDNDXINTERPOLANT_INSTANTIATION
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXInterpolant.h"
#include "PROPOSAL/methods.h"

#include <cmath>

#include "CubicInterpolation/Axis.h"
#include "CubicInterpolation/FindParameter.hpp"
using namespace PROPOSAL;

template <>
std::function<double(double, double, double)> transform_loss<crosssection::ComptonKleinNishina>::func
    = [](double v_cut, double v_max, double v) {
          return transform_loss_log(v_cut, v_max, v);
      };

template <>
std::function<double(double, double, double)> retransform_loss<crosssection::ComptonKleinNishina>::func
    = [](double v_cut, double v_max, double v) {
          return retransform_loss_log(v_cut, v_max, v);
      };

namespace PROPOSAL {
    double transform_relative_loss(double v_cut, double v_max, double v) {
        if (v < 0 || v_max == 0)
            return v_cut;
        if (v >= 1)
            return v_max;
        return v_cut * std::exp(v * std::log(v_max / v_cut));
    }

    double retransform_relative_loss(double v_cut, double v_max, double v) {
        if (v <= v_cut)
            return 0;
        if (v >= v_max)
            return 1;
        return std::log(v / v_cut) / std::log(v_max / v_cut);
    }

    double transform_loss_log(double v_cut, double v_max, double v)
    {
        auto xi = std::log((1. - v_cut)/(1 - v_max));
        return 1. - (1. - v_max) * std::exp((1. - v) * xi);
    }

    double retransform_loss_log(double v_cut, double v_max, double v)
    {
        auto xi = std::log((1. - v_cut)/(1 - v_max));
        return 1. - std::log((1. - v)/(1. - v_max)) / xi;
    }
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

size_t CrossSectionDNDXInterpolant::gen_hash(size_t hash) const {
    hash_combine(hash,
                 InterpolationSettings::NODES_DNDX_E,
                 InterpolationSettings::NODES_DNDX_V,
                 InterpolationSettings::UPPER_ENERGY_LIM);
    return hash;
}

double CrossSectionDNDXInterpolant::Calculate(double energy)
{
    auto lim = GetIntegrationLimits(energy);
    return Calculate(energy, lim.max);
}

double CrossSectionDNDXInterpolant::Calculate(double energy, double v)
{
    if (energy < lower_energy_lim)
        return 0.;
    auto lim = GetIntegrationLimits(energy);
    v = retransform_v(lim.min, lim.max, v);
    return interpolant.evaluate(std::array<double, 2> { energy, v });
}

double CrossSectionDNDXInterpolant::GetUpperLimit(double energy, double rate)
{
    if (energy < lower_energy_lim)
        throw std::invalid_argument("no dNdx for this energy defined.");
    auto lim = GetIntegrationLimits(energy);

    auto initial_guess = cubic_splines::ParameterGuess<std::array<double, 2>>();
    initial_guess.x = { energy, NAN };
    initial_guess.n = 1;

    auto v = cubic_splines::find_parameter(interpolant, rate, initial_guess);
    return transform_v(lim.min, lim.max, v);
}
