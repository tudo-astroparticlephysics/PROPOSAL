#include <algorithm>

#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/propagation_utility/Displacement.h"

using namespace PROPOSAL;
using std::max;

Displacement::Displacement(const CrossSectionList& cross)
    : cross(cross)
{
    if (cross.size() < 1)
        throw std::invalid_argument("at least one crosssection is required.");
}

double Displacement::FunctionToIntegral(
    const ParticleDef& p_def, const Medium& medium, double energy)
{
    auto result = 0.0;
    for (const auto& cr : cross)
        result += cr->CalculatedEdx(p_def, medium, energy);

    return -1.0 / result;
}

size_t Displacement::GetHash(const ParticleDef& p_def, const Medium& medium) const
{
    auto hash_digest = GetHash(p_def, medium);
    for (const auto& c : cross)
        hash_combine(hash_digest, c->GetHash(p_def, medium));
    return hash_digest;
}

double Displacement::GetLowerLim(const ParticleDef& p_def) const
{
    auto lower_lim = (double)0;
    for (const auto& c : cross)
        lower_lim = max(lower_lim, c->GetLowerEnergyLim(p_def));
    return lower_lim;
}

namespace PROPOSAL {
Interpolant1DBuilder::Definition displacement_interpol_def;
} // namespace PROPOSAL
