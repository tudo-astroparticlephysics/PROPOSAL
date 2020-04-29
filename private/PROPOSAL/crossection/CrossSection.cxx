#include <cassert>

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/methods.h"

using std::max;
using std::min;

using namespace PROPOSAL;

double CrossSection::GetEnergyCut(double energy) const
{
    auto limits = parametrization_->GetKinematicLimits(energy);

    if (cuts_)
        return min(max(limits.vMin, cuts_->GetCut(energy)), limits.vMax);
    else
        return limits.vMin;
}

size_t CrossSection::GetHash() const
{
    size_t hash_digest = 0;
    hash_combine(hash_digest, parametrization_->GetHash(), cuts_->GetHash());

    return hash_digest;
}

double CrossSection::GetLowerEnergyLimit() const
{
    return parametrization_->GetLowerEnergyLim();
};
