
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/propagation_utility/Displacement.h"

using namespace PROPOSAL;

Displacement::Displacement(const CrossSectionList& cross)
    : cross(cross)
    , lower_lim(std::numeric_limits<double>::max())
{
    if (cross.size() < 1)
        throw std::invalid_argument("at least one crosssection is required.");

    for (auto c : cross)
        lower_lim = std::min(lower_lim, c->GetParametrization().GetLowerEnergyLim());
}

double Displacement::FunctionToIntegral(double energy)
{
    auto result = 0.0;
    for (const auto& cr : cross)
        result += cr->CalculatedEdx(energy);

    return -1.0 / result;
}

namespace PROPOSAL {
Interpolant1DBuilder::Definition displacement_interpol_def(
    nullptr, 1000, 0., 1e14, 5, false, false, true, 5, false, false, false);
} // namespace PROPOSAL
