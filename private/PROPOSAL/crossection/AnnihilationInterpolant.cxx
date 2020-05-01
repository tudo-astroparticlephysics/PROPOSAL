#include "PROPOSAL/crossection/AnnihilationInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Annihilation.h"

using namespace PROPOSAL;

AnnihilationInterpolant::AnnihilationInterpolant(
    unique_ptr<Annihilation>&& param, const InterpolationDef& interpol_def)
    : CrossSectionInterpolant(forward<unique_ptr<Annihilation>>(param), nullptr,  interpol_def)
{
}

double AnnihilationInterpolant::CalculatedEdx(double energy)
{
    (void)energy;
    return 0;
}

double AnnihilationInterpolant::CalculatedE2dx(double energy)
{
    (void)energy;
    return 0;
}
