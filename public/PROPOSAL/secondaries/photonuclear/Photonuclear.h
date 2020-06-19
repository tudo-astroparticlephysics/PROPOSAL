#pragma once
#include <stdexcept>

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/DefaultFactory.h"

namespace PROPOSAL {
namespace secondaries {
    struct Photonuclear : public Parametrization,
                          public RegisteredInDefault<Photonuclear> {
        Photonuclear(ParticleDef, Medium){};

        static constexpr size_t n_rnd = 0;
        size_t RequiredRandomNumbers() { return n_rnd; }

        static constexpr InteractionType type
            = PROPOSAL::InteractionType::Photonuclear;
        InteractionType GetInteractionType() { return type; }

        vector<Loss::secondary_t> CalculateSecondaries(
            double, Loss::secondary_t, const Component&, vector<double>)
        {
            throw std::logic_error(
                "decay objects of a photonuclear interaction cannot be "
                "calculated by PROPOSAL.");
        };
    };
} // namespace secondaries
} // namespace PROPOSAL
