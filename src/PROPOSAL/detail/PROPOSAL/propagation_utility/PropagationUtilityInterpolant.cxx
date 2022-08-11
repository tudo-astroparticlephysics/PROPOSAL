
#include <cassert>
#include <cmath>
#include <functional>

#include "CubicInterpolation/Interpolant.h"
#include "CubicInterpolation/FindParameter.hpp"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/MathMethods.h"

using namespace PROPOSAL;

std::string UtilityInterpolant::gen_path() const
{
    return std::string(InterpolationSettings::TABLES_PATH);
}

std::string UtilityInterpolant::gen_name(std::string prefix) const
{
    return std::string(prefix) + std::to_string(this->GetHash())
        + std::string(".dat");
}
UtilityInterpolant::UtilityInterpolant(
    std::function<double(double)> f, double lim, size_t hash)
    : UtilityIntegral(f, lim, hash)
    , lower_lim(lim)
    , interpolant_(nullptr)
{
}

void UtilityInterpolant::BuildTables(const std::string prefix, size_t nodes,
                                     bool reverse) {
    auto def = cubic_splines::CubicSplines<double>::Definition();
    auto reference_x = lower_lim;
    reverse_ = reverse;
    if (reverse_) {
        reference_x = InterpolationSettings::UPPER_ENERGY_LIM;
    }

    hash_combine(this->hash, nodes, reverse,
                 InterpolationSettings::UPPER_ENERGY_LIM);

    if (reverse) {
        def.f = [&](double energy) {
            return UtilityIntegral::Calculate(reference_x, energy);
        };
    } else {
        def.f = [&](double energy) {
            return UtilityIntegral::Calculate(energy, reference_x);
        };
    }
    def.f_trafo = std::make_unique<cubic_splines::ExpM1Axis<double>>(1., 0.);
    def.axis = std::make_unique<cubic_splines::ExpAxis<double>>(
            lower_lim, InterpolationSettings::UPPER_ENERGY_LIM, nodes);

    interpolant_ = std::make_shared<interpolant_t>(
            std::move(def), gen_path(), gen_name(prefix));
}


double UtilityInterpolant::Calculate(double energy_initial, double energy_final)
{
    assert(energy_initial >= energy_final);
    assert(energy_final >= lower_lim);

    if (energy_initial - energy_final < energy_initial * IPREC)
        return FunctionToIntegral((energy_initial + energy_initial) / 2)
            * (energy_final - energy_initial);

    auto integral_upper_limit = interpolant_->evaluate(energy_initial);
    auto integral_lower_limit = interpolant_->evaluate(energy_final);

    if (reverse_)
        return integral_lower_limit - integral_upper_limit;
    return integral_upper_limit - integral_lower_limit;
}

// ------------------------------------------------------------------------- //
double UtilityInterpolant::GetUpperLimit(double upper_limit, double rnd)
{
    assert(rnd >= 0);

    auto max_rnd = Calculate(upper_limit, lower_lim);
    if (rnd > max_rnd)
        throw std::logic_error("Unable to calculate GetUpperLimit since result"
                               "is below lower_lim. rnd was " + std::to_string(rnd)
                               + " with rnd_max " + std::to_string(max_rnd));

    if (rnd == max_rnd)
        return lower_lim;

    if (reverse_)
        rnd = -rnd;

    auto integrated_to_upper = interpolant_->evaluate(upper_limit);
    auto initial_guess = cubic_splines::ParameterGuess<double>();

    // find initial parameters for newton raphson method by using bisection
    auto f = [this, &integrated_to_upper, &rnd](double val) {
        return interpolant_->evaluate(val) - (integrated_to_upper - rnd);
    };
    auto bisec_tolerance = (upper_limit - lower_lim) * 1e-2;
    std::tie(initial_guess.lower, initial_guess.upper) =
            Bisection(f, lower_lim, upper_limit, bisec_tolerance, 100);

    // if we are close to the upper limit, start newton raphson method there
    if (initial_guess.upper == upper_limit) {
        initial_guess.x = upper_limit;
    }
    else
        initial_guess.x = (initial_guess.lower + initial_guess.upper) / 2;

    try {
        return cubic_splines::find_parameter(
                *interpolant_, integrated_to_upper - rnd, initial_guess);
    } catch (std::runtime_error&) {
        Logging::Get("proposal.UtilityInterpolant")->warn(
                "Newton-Raphson iteration in UtilityInterpolant::GetUpperLimit "
                "failed. Try solving using bisection method.");

        return Bisection(f, lower_lim, upper_limit, 1e-6, 100).first;
    }

    // TODO: Check whether this is already accurate enough
    // (see e.g. version at a81e54f62f4383936cb046da4cad7429a48bb750 for old
    // version)
}
