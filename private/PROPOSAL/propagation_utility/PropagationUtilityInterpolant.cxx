
#include <cassert>
#include <cmath>
#include <functional>

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/methods.h"

#include "PROPOSAL/crossection/CrossSection.h"

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
    Interpolant1DBuilder::Definition interpol_def)
{
    auto utility_func = [&](double energy) {
        return UtilityIntegral::Calculate(energy, lower_lim);
    };
    interpol_def.function1d = utility_func;
    interpol_def.xmin = lower_lim;
    auto builder = make_unique<Interpolant1DBuilder>(interpol_def);
    interpolant_ = InitializeInterpolation(
        name, std::move(builder), hash, InterpolationDef());
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

    assert(integrated_to_upper > rnd); // searched Energy is below lower_lim
                                       // return lower_lim as a lower limit

    if (upper_limit - lower_limit > upper_limit * IPREC)
        return lower_limit;

    auto step = upper_limit + 0.5 * rnd / FunctionToIntegral(upper_limit);

    return upper_limit + rnd / FunctionToIntegral(step);
}
