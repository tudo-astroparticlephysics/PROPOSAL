#include <algorithm>

#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/propagation_utility/Displacement.h"

using namespace PROPOSAL;
using std::max;

template <typename Cross>
double Displacement::FunctionToIntegral(Cross&& cross, double energy)
{
    auto result = 0.0;
    for (auto& cr : cross)
        result += cr->CalculatedEdx(energy);

    return -1.0 / result;
}

template <typename Cross> size_t Displacement::GetHash(Cross&& cross) const
{
    auto hash_digest = size_t{ 0 };
    for (const auto& c : cross)
        hash_combine(hash_digest, c->GetHash());
    return hash_digest;
}
template <typename Cross> double Displacement::GetLowerLim(Cross&& cross) const
{
    auto lower_lim = (double)0;
    for (const auto& c : cross)
        lower_lim = max(lower_lim, c->GetLowerEnergyLim());
    return lower_lim;
}

namespace PROPOSAL {
Interpolant1DBuilder::Definition displacement_interpol_def;
} // namespace PROPOSAL
