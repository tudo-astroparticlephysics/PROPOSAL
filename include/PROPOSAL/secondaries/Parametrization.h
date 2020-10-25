#pragma once

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/particle/Particle.h"
#include <vector>

#include <iostream>

using PROPOSAL::Components::Component;
using std::vector;

namespace PROPOSAL {
namespace secondaries {
    struct Parametrization {
        Parametrization() = default;
        virtual ~Parametrization() = default;

        virtual size_t RequiredRandomNumbers() const noexcept = 0;
        virtual InteractionType GetInteractionType() const noexcept = 0;
        virtual vector<DynamicData> CalculateSecondaries(StochasticLoss,
                                    const Component&, vector<double> &rnd) = 0;
    };
} // namespace secondaries
} // namespace PROPOSAL
