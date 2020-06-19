#pragma once

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/medium/Components.h"
#include <vector>

#include <iostream>

using std::vector;
using PROPOSAL::Components::Component;

namespace PROPOSAL {
namespace secondaries {
    struct Parametrization {
        Parametrization() = default;
        ~Parametrization() = default;

        virtual size_t RequiredRandomNumbers() = 0;
        virtual vector<Loss::secondary_t> CalculateSecondaries(
            double primary_energy, Loss::secondary_t, const Component&,
            vector<double> rnd)
            = 0;
        virtual InteractionType GetInteractionType() = 0;
    };
} // namespace secondaries
} // namespace PROPOSAL
