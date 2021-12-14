#pragma once

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/particle/ParticleDef.h"

namespace PROPOSAL {
namespace secondaries {
    struct Photoproduction : public Parametrization,
            public DefaultSecondaries<Photoproduction> {
        public:
        Photoproduction(const ParticleDef&, const Medium&) {};

        static constexpr size_t n_rnd = 0;
        size_t RequiredRandomNumbers() const noexcept override { return n_rnd; }

        static constexpr InteractionType type
        = PROPOSAL::InteractionType::Photoproduction;
        InteractionType GetInteractionType() const noexcept override { return type; }

        std::vector<ParticleState> CalculateSecondaries(
                StochasticLoss loss, const Component&, std::vector<double>&) override
                {
            // all energies goes to hadronic component
            ParticleState hadron;
            hadron.energy = loss.energy;
            hadron.SetType(ParticleType::Hadron);
            hadron.time = loss.time;
            hadron.position = loss.position;
            hadron.direction = loss.direction;
            hadron.propagated_distance = 0.;

            return std::vector<ParticleState>{hadron};
                };
    };
} // namespace secondaries
} // namespace PROPOSAL
