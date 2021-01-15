#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDX.h"

using namespace PROPOSAL;

std::array<double, 2> CrossSectionDNDX::GetIntegrationLimits(
    double energy) const
{
    auto kin_lim = kinematic_limits(energy);
    auto lim = std::array<double, 2> {
        std::get<crosssection::Parametrization::V_MIN>(kin_lim),
        std::get<crosssection::Parametrization::V_MAX>(kin_lim)
    };
    if (cut)
        lim[0] = std::max(cut->GetCut(kin_lim, energy), lim[0]);
    return lim;
}
