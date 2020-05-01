
#include "PROPOSAL/crossection/AnnihilationIntegral.h"
#include "PROPOSAL/crossection/parametrization/Annihilation.h"

using namespace PROPOSAL;

AnnihilationIntegral::AnnihilationIntegral(unique_ptr<Annihilation>&& param)
    : CrossSectionIntegral(forward<unique_ptr<Annihilation>>(param), nullptr)
{
}

double AnnihilationIntegral::CalculatedEdx(double energy)
{
    (void)energy;
    return 0;
}

double AnnihilationIntegral::CalculatedE2dx(double energy)
{
    (void)energy;
    return 0;
}
