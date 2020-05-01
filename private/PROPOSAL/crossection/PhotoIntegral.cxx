#include "PROPOSAL/crossection/PhotoIntegral.h"
#include "PROPOSAL/crossection/parametrization/Photonuclear.h"

using namespace PROPOSAL;

PhotoIntegral::PhotoIntegral(unique_ptr<Photonuclear>&& param, shared_ptr<const EnergyCutSettings> cut)
    : CrossSectionIntegral(forward<unique_ptr<Photonuclear>>(param), cut) {}
