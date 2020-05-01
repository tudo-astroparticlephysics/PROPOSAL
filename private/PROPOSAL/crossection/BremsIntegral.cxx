
#include <functional>

#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crossection/BremsIntegral.h"

using namespace PROPOSAL;

BremsIntegral::BremsIntegral(unique_ptr<Bremsstrahlung>&& param, shared_ptr<const EnergyCutSettings> cut)
    : CrossSectionIntegral(forward<unique_ptr<Bremsstrahlung>>(param), cut)
{}
