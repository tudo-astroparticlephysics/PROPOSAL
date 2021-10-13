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

        std::vector<ParticleState> CalculateSecondaries(
                StochasticLoss loss, const Component&, std::vector<double>&) override
        {
            ParticleState primary_lepton;
            primary_lepton.energy = loss.parent_particle_energy - loss.energy;
            primary_lepton.type = primary_particle_type;
            primary_lepton.time = loss.time;
            primary_lepton.position = loss.position;
            primary_lepton.direction = loss.direction;
            primary_lepton.propagated_distance = 0.;

            ParticleState hadron;
            hadron.energy = loss.energy;
            hadron.SetType(ParticleType::Hadron);
            hadron.time = loss.time;
            hadron.position = loss.position;
            hadron.direction = loss.direction;
            hadron.propagated_distance = 0.;

            return std::vector<ParticleState>{primary_lepton, hadron};
        };
    };
} // namespace secondaries
} // namespace PROPOSAL
