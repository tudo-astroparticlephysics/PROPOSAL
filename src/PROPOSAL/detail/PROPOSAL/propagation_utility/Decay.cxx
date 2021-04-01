#include "PROPOSAL/propagation_utility/Decay.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/propagation_utility/Displacement.h"

using namespace PROPOSAL;

double Decay::FunctionToIntegral(double energy) {
    assert(!std::isinf(lifetime));
    assert(energy >= mass);
    double square_momentum = (energy - mass) * (energy + mass);
    double aux = SPEED * std::sqrt(square_momentum) / mass;
    return disp->FunctionToIntegral(energy) / aux;
}

Decay::Decay(std::shared_ptr<Displacement> _disp, double _lifetime, double _mass)
    : lifetime(_lifetime)
            , mass(_mass)
            , disp(_disp)
            , hash(0)
    {
        hash_combine(hash, disp->GetHash(), lifetime, mass);
    }
