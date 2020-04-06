#include "PROPOSAL/propagation_utility/ContRand.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

using namespace PROPOSAL;

ContRand::ContRand(CrossSectionList cross)
    : cross(cross)
    , mass(cross.front()->GetParametrization().GetParticleMass())
{
}

namespace PROPOSAL {
Interpolant1DBuilder::Definition contrand_interpol_def(
    nullptr, 200, 0., 1e14, 5, false, false, true, 5, false, false, false);
} // namespace PROPOSAL
