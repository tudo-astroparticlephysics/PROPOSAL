#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDX.h"

using namespace PROPOSAL;

std::array<double, 2> CrossSectionDNDX::GetIntegrationLimits(
    double energy) const
{
    auto kin_lim = kinematic_limits(energy);
    auto lim = std::array<double, 2>();
    if (cut)
        lim[0] = cut->GetCut(kin_lim, energy);
    else
        lim[0] = std::get<crosssection::Parametrization::V_MIN>(kin_lim);
    lim[1] = std::get<crosssection::Parametrization::V_MAX>(kin_lim);
    return lim;
}

double CrossSectionDNDX::GetLowerEnergyLim() const { return lower_energy_lim; }
