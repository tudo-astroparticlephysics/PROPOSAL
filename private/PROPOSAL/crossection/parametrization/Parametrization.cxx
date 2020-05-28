
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include <sstream>

using namespace PROPOSAL;
using std::string;

Parametrization::Parametrization(InteractionType type, const string& name)
    : interaction_type(type)
    , name(name)
{
}

Parametrization::KinematicLimits::KinematicLimits(double v_min, double v_max)
    : vMin(v_min)
    , vMax(v_max)
{
}

size_t Parametrization::GetHash() const
{
    size_t hash_digest = 0;
    hash_combine(hash_digest, name);

    return hash_digest;
}
