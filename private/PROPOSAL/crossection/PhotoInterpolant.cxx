#include <functional>

#include "PROPOSAL/crossection/PhotoInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Photonuclear.h"

using namespace PROPOSAL;

PhotoInterpolant::PhotoInterpolant(unique_ptr<Photonuclear>&& param,
    shared_ptr<const EnergyCutSettings> cut, const InterpolationDef& def)
    : CrossSectionInterpolant(forward<unique_ptr<Photonuclear>>(param), cut, def)
{
}
