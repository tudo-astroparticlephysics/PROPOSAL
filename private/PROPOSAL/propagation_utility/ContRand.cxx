#include "PROPOSAL/propagation_utility/ContRand.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

using namespace PROPOSAL;

ContRand::ContRand(CrossSectionList cross)
    : cross(cross)
    , lower_lim(std::numeric_limits<double>::max())
{
    if (cross.size() < 1)
        throw std::invalid_argument("at least one crosssection is required.");

    for (auto c : cross)
        lower_lim = std::min(lower_lim, c->GetParametrization().GetLowerEnergyLim());
}

namespace PROPOSAL {
Interpolant1DBuilder::Definition contrand_interpol_def(
    nullptr, 200, 0., 1e14, 5, false, false, true, 5, false, false, false);
} // namespace PROPOSAL
