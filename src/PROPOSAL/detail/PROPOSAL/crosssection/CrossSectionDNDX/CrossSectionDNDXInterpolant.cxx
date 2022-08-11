#define CROSSSECTIONDNDXINTERPOLANT_INSTANTIATION
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXInterpolant.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>

#include "CubicInterpolation/Axis.h"
#include "CubicInterpolation/FindParameter.hpp"
using namespace PROPOSAL;

namespace PROPOSAL {
    double transform_relative_loss(double v_cut, double v_max, double v, double c) {
        if (v < 0 || v_max == 0)
            return v_cut;
        if (v >= 1)
            return v_max;
        return v_cut * std::exp(std::pow(v, c) * std::log(v_max / v_cut));
    }

    double retransform_relative_loss(double v_cut, double v_max, double v, double c) {
        if (v <= v_cut)
            return 0;
        if (v >= v_max)
            return 1;
        return std::pow(std::log(v / v_cut) / std::log(v_max / v_cut), 1. / c);
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
        + std::string(".dat");
}

size_t CrossSectionDNDXInterpolant::gen_hash(size_t hash) const {
    hash_combine(hash,
                 InterpolationSettings::NODES_DNDX_E,
                 InterpolationSettings::NODES_DNDX_V,
                 InterpolationSettings::UPPER_ENERGY_LIM);
    return hash;
}

double CrossSectionDNDXInterpolant::evaluate_interpolant(double E, double vbar)
{
    if (E < lower_energy_lim)
        return 0.;
    auto dNdx = interpolant.evaluate(std::array<double, 2> { E, vbar });
    if (dNdx < 0) {
        auto inter_name = Type_Interaction_Name_Map.at(type_id);
        logger->warn("Negative dNdx value for E = {:.4f} MeV, vbar = {:.4f} "
                     "detected in interaction type {}. "
                     "Setting dNdx to zero.", E, vbar, inter_name);
        dNdx = 0.;
    }
    return dNdx;
}

double CrossSectionDNDXInterpolant::Calculate(double energy)
{
    // for v=v_max, we know by construction of the transformation that vbar=1
    return evaluate_interpolant(energy, 1);
}

double CrossSectionDNDXInterpolant::Calculate(double energy, double v)
{
    auto lim = GetIntegrationLimits(energy);
    auto vbar = retransform_v(lim.min, lim.max, v);
    return evaluate_interpolant(energy, vbar);
}

double CrossSectionDNDXInterpolant::GetUpperLimit(double energy, double rate)
{
    if (energy < lower_energy_lim)
        throw std::invalid_argument("no dNdx for this energy defined.");
    auto lim = GetIntegrationLimits(energy);

    auto initial_guess = cubic_splines::ParameterGuess<std::array<double, 2>>();
    initial_guess.x = { energy, NAN };
    initial_guess.n = 1;
    double v;
    try {
        v = cubic_splines::find_parameter(interpolant, rate, initial_guess);
    } catch (std::runtime_error&) {
        Logging::Get("proposal.UtilityInterpolant")->warn(
                "Newton-Raphson iteration in "
                "CrossSectionDNDXInterpolant::GetUpperLimit failed. Try solving"
                " using bisection method.");

        auto f = [this, &rate, &energy](double val) {
            return interpolant.evaluate(
                    std::array<double, 2> { energy, val }) - rate;
        };
        // v is evaluated in transformed space!
        auto interval = Bisection(f, retransform_v(lim.min, lim.max, lim.min),
                                  retransform_v(lim.min, lim.max, lim.max),
                                  1e-6, 100);
        v = (interval.first + interval.second) / 2.;
    }
    return transform_v(lim.min, lim.max, v);
}
