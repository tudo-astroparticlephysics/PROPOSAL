
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include <sstream>

using namespace PROPOSAL;
using std::make_tuple;
using std::string;

Parametrization::Parametrization(
    InteractionType interaction_type, const string& name)
    : interaction_type(interaction_type)
    , name(name)
{
}

double Parametrization::DifferentialCrossSection(
    const ParticleDef&, const Component&, double, double)
{
    throw std::logic_error("Not implemented error.");
}

tuple<double, double> Parametrization::GetKinematicLimits(
    const ParticleDef&, const Component&, double) const noexcept
{
    throw std::logic_error("Not implemented error.");
}

size_t Parametrization::GetHash() const noexcept
{
    auto hash_digest = size_t{ 0 };
    hash_combine(hash_digest, name);
    return hash_digest;
}
