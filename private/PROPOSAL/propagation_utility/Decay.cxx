
#include "PROPOSAL/propagation_utility/Decay.h"

using namespace PROPOSAL;

Decay::Decay(const CrossSectionList& cross)
    : cross(cross)
    , mass(cross.front()->GetParametrization().GetParticleMass())
    , lifetime(cross.front()->GetParametrization().GetParticleLifetime())
{
}

namespace PROPOSAL {
Interpolant1DBuilder::Definition decay_interpol_def(
    nullptr, 200, 0., 1e14, 5, false, false, true, 5, false, false, false);
} // namespace PROPOSAL
