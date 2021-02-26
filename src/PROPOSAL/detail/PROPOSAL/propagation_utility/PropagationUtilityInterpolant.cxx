
#include <cassert>
#include <cmath>
#include <functional>

#include "CubicInterpolation/Interpolant.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

using namespace PROPOSAL;
using Helper::InitializeInterpolation;

std::string UtilityInterpolant::gen_path() const
{
    return std::string(InterpolationSettings::TABLES_PATH);
}

std::string UtilityInterpolant::gen_name(std::string prefix) const
{
    return std::string(prefix) + std::to_string(this->GetHash())
        + std::string(".txt");
}
UtilityInterpolant::UtilityInterpolant(
    std::function<double(double)> f, double lim, size_t hash)
    : UtilityIntegral(f, lim, hash)
    , lower_lim(lim)
    , interpolant_(nullptr)
{
}

/* void UtilityInterpolant::BuildTables(std::string name, */
/*     Interpolant1DBuilder::Definition interpol_def, bool reverse) */
/* { */
/* interpol_def.xmin = lower_lim; */
/* double reference_x; */
/* if (reverse) { */
/*     reference_x = interpol_def.xmax; */
/* } else { */
/*     reference_x = interpol_def.xmin; */
/* } */
/* } */

double UtilityInterpolant::Calculate(double energy_initial, double energy_final)
{
    assert(energy_initial >= energy_final);
    assert(energy_final >= lower_lim);

    if (energy_initial - energy_final < energy_initial * IPREC)
        return FunctionToIntegral((energy_initial + energy_initial) / 2)
            * (energy_final - energy_initial);

    auto integral_upper_limit = interpolant_->evaluate(energy_initial);
    auto integral_lower_limit = interpolant_->evaluate(energy_final);

    return integral_upper_limit - integral_lower_limit;
}

// ------------------------------------------------------------------------- //
double UtilityInterpolant::GetUpperLimit(double upper_limit, double rnd)
{
    assert(rnd >= 0);

    auto integrated_to_upper = interpolant_->evaluate(upper_limit);
    auto x_guess = cubic_splines::find_parameter(
        *interpolant_, integrated_to_upper - rnd, rnd);

    assert(integrated_to_upper >= rnd
        or integrated_to_upper < 0); // searched Energy is below lower_lim
                                     // or we use reverse interpolation

    if (std::abs(upper_limit - lower_lim) > upper_limit * IPREC)
        return lower_lim;

    return x_guess;
    // TODO: Check whether this is already accurate enough
    /*
    auto step = upper_limit + 0.5 * rnd / FunctionToIntegral(upper_limit);

    return upper_limit + rnd / FunctionToIntegral(step);
     */
}
