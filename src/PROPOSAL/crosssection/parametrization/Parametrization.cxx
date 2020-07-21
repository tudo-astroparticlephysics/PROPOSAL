
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include <sstream>

using namespace PROPOSAL;
using std::make_tuple;
using std::string;

crosssection::Parametrization::Parametrization(
    InteractionType interaction_type, const string& name)
    : interaction_type(interaction_type)
    , name(name)
{
}

double crosssection::Parametrization::DifferentialCrossSection(
    const ParticleDef&, const Component&, double, double) const
{
    throw std::logic_error("Not implemented error.");
}

tuple<double, double> crosssection::Parametrization::GetKinematicLimits(
    const ParticleDef&, const Component&, double) const noexcept
{
    throw std::logic_error("Not implemented error.");
}

size_t crosssection::Parametrization::GetHash() const noexcept
{
    auto hash_digest = size_t{ 0 };
    hash_combine(hash_digest, name);
    return hash_digest;
}
