
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
{
}

void UtilityInterpolant::BuildTables(std::string name, size_t hash, Interpolant1DBuilder::Definition interpol_def)
{
    Interpolant1DBuilder interpolant_builder(interpol_def);
    interpolant_ = InitializeInterpolation(name, interpolant_builder, hash, InterpolationDef());
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
