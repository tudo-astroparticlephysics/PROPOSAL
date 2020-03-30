
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
/******************************************************************************
 *                              Utility Integral                              *
 ******************************************************************************/

UtilityInterpolant::UtilityInterpolant(std::function<double(double)> func)
    : UtilityIntegral(func)
    , interpolant_(nullptr)
    /* , interpolant_diff_(nullptr) */
    , low_(1e-3)
{
    builder1d.SetMax(utility_interpolation_def.nodes_propagate)
        .SetXMin(low_)
        .SetXMax(utility_interpolation_def.max_node_energy)
        .SetRomberg(utility_interpolation_def.order_of_interpolation)
        .SetRational(false)
        .SetRelative(false)
        .SetIsLog(true)
        .SetRombergY(utility_interpolation_def.order_of_interpolation)
        .SetRationalY(false)
        .SetRelativeY(false)
        .SetLogSubst(false)
        .SetFunction1D([this](double energy){ return UtilityIntegral::Calculate(energy, low_, 0);});

    /* builder_diff.SetMax(utility_interpolation_def.nodes_propagate) */
    /*     .SetXMin(low_) */
    /*     .SetXMax(utility_interpolation_def.max_node_energy) */
    /*     .SetRomberg(utility_interpolation_def.order_of_interpolation) */
    /*     .SetRational(false) */
    /*     .SetRelative(false) */
    /*     .SetIsLog(true) */
    /*     .SetRombergY(utility_interpolation_def.order_of_interpolation) */
    /*     .SetRationalY(false) */
    /*     .SetRelativeY(false) */
    /*     .SetLogSubst(false) */
    /*     .SetFunction1D(FunctionToIntegral); */
}

void UtilityInterpolant::BuildTables(std::string name, size_t hash)
{
    interpolant_ = InitializeInterpolation(name, builder1d, hash, utility_interpolation_def);

    /* interpolant_diff_ = InitializeInterpolation( */
    /*     name.append("diff"), builder_diff, hash, utility_interpolation_def); */
}

double UtilityInterpolant::Calculate(double ei, double ef, double rnd)
{
    (void)rnd;
    assert(ef > 0);
    assert(ei > ef);

    upper_limit = std::make_pair(interpolant_->Interpolate(ei), ei);

    if (ei - ef < ei * IPREC)
        return upper_limit.first * (ei - ef);

    return upper_limit.first - interpolant_->Interpolate(ef);
}

// ------------------------------------------------------------------------- //
double UtilityInterpolant::GetUpperLimit(double ei, double rnd)
{
    if (ei != upper_limit.second)
        Calculate(ei, 1e-3, rnd);

    auto lower_limit = interpolant_->FindLimit(upper_limit.first - rnd);

    if (upper_limit.first - lower_limit > upper_limit.first * IPREC)
        return lower_limit;

    auto initial_step = ei + 0.5 * rnd / FunctionToIntegral(ei);
    return ei + rnd / FunctionToIntegral(initial_step);
}

InterpolationDef UtilityInterpolant::utility_interpolation_def;
