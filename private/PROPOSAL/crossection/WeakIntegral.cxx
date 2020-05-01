#include "PROPOSAL/crossection/WeakIntegral.h"
#include "PROPOSAL/crossection/parametrization/WeakInteraction.h"

using namespace PROPOSAL;

WeakIntegral::WeakIntegral(unique_ptr<WeakInteraction>&& param)
    : CrossSectionIntegral(forward<unique_ptr<WeakInteraction>>(param), nullptr)
{
}

double WeakIntegral::CalculatedEdx(double energy)
{
    (void)energy;
    return 0;
}

double WeakIntegral::CalculatedE2dx(double energy)
{
    (void)energy;
    return 0;
}
