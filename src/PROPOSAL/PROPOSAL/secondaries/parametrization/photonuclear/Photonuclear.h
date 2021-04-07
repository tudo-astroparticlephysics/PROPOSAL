#pragma once
#include <stdexcept>

#include "PROPOSAL/particle/Particle.h"

namespace PROPOSAL {
namespace secondaries {
    struct Photonuclear : public Parametrization,
                          public DefaultSecondaries<Photonuclear> {
        int primary_particle_type;

        Photonuclear(ParticleDef p, Medium)
            : primary_particle_type(p.particle_type){};

        static constexpr size_t n_rnd = 0;
        size_t RequiredRandomNumbers() const noexcept override { return n_rnd; }

        static constexpr InteractionType type
            = PROPOSAL::InteractionType::Photonuclear;
        InteractionType GetInteractionType() const noexcept override { return type; }

        std::vector<ParticleState> CalculateSecondaries(StochasticLoss, const Component&,
                                                   std::vector<double>&) override
        {
            Logging::Get("proposal.Secondaries")
                    ->warn("PROPOSAL can not generate secondary particles"
                            "for a photonuclear interaction.");
            //TODO: Treatment for hadronic interactions
            return std::vector<ParticleState>{};
        };
    };
} // namespace secondaries
} // namespace PROPOSAL
