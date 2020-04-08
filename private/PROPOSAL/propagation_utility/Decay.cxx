
#include "PROPOSAL/propagation_utility/Decay.h"

using namespace PROPOSAL;

Decay::Decay(const CrossSectionList& cross, const ParticleDef& p_def) : Decay(cross, p_def.lifetime, p_def.mass){};

Decay::Decay(const CrossSectionList& cross, double lifetime, double mass)
    : cross(cross)
    , mass(mass)
    , lifetime(lifetime)
    , lower_lim(std::numeric_limits<double>::max())
{
    if (cross.size() < 1)
        throw std::invalid_argument("at least one crosssection is required.");

    for (auto c : cross)
        lower_lim = std::min(lower_lim, c->GetParametrization().GetLowerEnergyLim());
}

namespace PROPOSAL {
Interpolant1DBuilder::Definition decay_interpol_def(
    nullptr, 200, 0., 1e14, 5, false, false, true, 5, false, false, false);
} // namespace PROPOSAL
