
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
    , interpolant_diff_(nullptr)
    , stored_result_(0)
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
        .SetFunction1D(std::bind(
            &UtilityIntegral::Calculate, this, std::placeholders::_1, low_, 0));

    builder_diff.SetMax(utility_interpolation_def.nodes_propagate)
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
        .SetFunction1D(FunctionToIntegral);
}

void UtilityInterpolant::BuildTables(std::string name, size_t hash)
{
    interpolant_ = InitializeInterpolation(name, builder1d, hash, utility_interpolation_def);

    interpolant_diff_ = InitializeInterpolation(name.append("diff"), builder_diff, hash, utility_interpolation_def);
}

// ------------------------------------------------------------------------- //
double UtilityInterpolant::GetUpperLimit(double ei, double rnd)
{
    return std::min(
        std::max(ei
                + rnd
                    / interpolant_diff_->Interpolate(
                          ei + rnd / (2 * interpolant_diff_->Interpolate(ei))),
            low_),
        ei);
}
