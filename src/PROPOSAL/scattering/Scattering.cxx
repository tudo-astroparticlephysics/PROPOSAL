#include "PROPOSAL/scattering/Scattering.h"

namespace PROPOSAL {
std::vector<std::unique_ptr<stochastic_deflection::Parametrization>>
make_stochastic_deflection(std::vector<InteractionType> const& types,
    ParticleDef const& p, Medium const& m)
{
    auto deflections = std::vector<
        std::unique_ptr<stochastic_deflection::Parametrization>>();
    for (auto t : types)
        deflections.emplace_back(
            DefaultFactory<stochastic_deflection::Parametrization>::Create(
                t, p, m));
    return deflections;
}
}// namespace PROPOSAL
