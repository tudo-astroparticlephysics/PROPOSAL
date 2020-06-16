#include "PROPOSAL/crossection/CrossSectionDNDX/CrossSectionDNDX.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"

using std::make_tuple;

using namespace PROPOSAL;

tuple<double, double> CrossSectionDNDX::GetIntegralLimits(
    EnergyCutSettings* cut, double energy)
{
    auto lim = kinematic_limits(energy);
    return make_tuple(
        cut->GetCut(lim, energy), get<Parametrization::V_MAX>(lim));
}

tuple<double, double> CrossSectionDNDX::GetIntegralLimits(
    std::nullptr_t, double energy)
{
    return kinematic_limits(energy);
}
