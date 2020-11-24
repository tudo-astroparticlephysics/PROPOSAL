
#include <cassert>
#include <cmath>
#include <functional>

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/methods.h"

#include "PROPOSAL/crosssection/CrossSection.h"

#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/math/MathMethods.h"

using namespace PROPOSAL;
using Helper::InitializeInterpolation;

UtilityInterpolant::UtilityInterpolant(
    std::function<double(double)> func, double lower_lim)
    : UtilityIntegral(func, lower_lim)
    , interpolant_(nullptr)
{
}

void UtilityInterpolant::BuildTables(std::string name, size_t hash,
    Interpolant1DBuilder::Definition interpol_def, bool reverse)
{
    interpol_def.xmin = lower_lim;
    double reference_x;
    if (reverse) {
        reference_x = interpol_def.xmax;
    } else {
        reference_x = interpol_def.xmin;
    }
    auto utility_func = [&](double energy) {
        return UtilityIntegral::Calculate(energy, reference_x);
    };
    interpol_def.function1d = utility_func;
    hash_combine(hash, reverse);
    interpolant_ = InitializeInterpolation(name, Interpolant1DBuilder(interpol_def), hash);
}

double UtilityInterpolant::Calculate(double energy_initial, double energy_final)
{
    assert(energy_initial >= energy_final);
    assert(energy_final >= lower_lim);

    if (energy_initial - energy_final < energy_initial * IPREC)
        return FunctionToIntegral((energy_initial + energy_initial) / 2)
            * (energy_final - energy_initial);

    auto integral_upper_limit = interpolant_->Interpolate(energy_initial);
    auto integral_lower_limit = interpolant_->Interpolate(energy_final);

    return integral_upper_limit - integral_lower_limit;
}

// ------------------------------------------------------------------------- //
double UtilityInterpolant::GetUpperLimit(double upper_limit, double rnd)
{
    assert(rnd >= 0);

    auto integrated_to_upper = interpolant_->Interpolate(upper_limit);
    auto lower_limit = interpolant_->FindLimit(integrated_to_upper - rnd);

    assert(integrated_to_upper >= rnd or integrated_to_upper < 0); // searched Energy is below lower_lim
                                                                  // or we use reverse interpolation

    if (std::abs(upper_limit - lower_limit) > upper_limit * IPREC)
        return lower_limit;

    auto step = upper_limit + 0.5 * rnd / FunctionToIntegral(upper_limit);

    return upper_limit + rnd / FunctionToIntegral(step);
}
