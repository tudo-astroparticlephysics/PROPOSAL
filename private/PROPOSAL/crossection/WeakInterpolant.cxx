
#include <functional>

#include "PROPOSAL/crossection/WeakInterpolant.h"
#include "PROPOSAL/crossection/parametrization/WeakInteraction.h"

using namespace PROPOSAL;

WeakInterpolant::WeakInterpolant(
    unique_ptr<WeakInteraction>&& param, const InterpolationDef& interpol_def)
    : CrossSectionInterpolant(forward<unique_ptr<WeakInteraction>>(param), nullptr, interpol_def)
{
}
double WeakInterpolant::CalculatedEdx(double energy)
{
    (void)energy;
    return 0;
}

double WeakInterpolant::CalculatedE2dx(double energy)
{
    (void)energy;
    return 0;
}
