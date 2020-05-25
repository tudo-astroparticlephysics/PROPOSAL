
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

double Parametrization::FunctionToDNdxIntegral(
    const ParticleDef& p_def, const Component& comp, double energy, double v)
{
    return DifferentialCrossSection(p_def, comp, energy, v);
}

double Parametrization::FunctionToDEdxIntegral(
    const ParticleDef& p_def, const Component& comp, double energy, double v)
{
    return v * DifferentialCrossSection(p_def, comp, energy, v);
}

double Parametrization::FunctionToDE2dxIntegral(
    const ParticleDef& p_def, const Component& comp, double energy, double v)
{
    return v * v * DifferentialCrossSection(p_def, comp, energy, v);
}

size_t Parametrization::GetHash() const
{
    size_t hash_digest = 0;
    hash_combine(hash_digest, name);

    return hash_digest;
}
