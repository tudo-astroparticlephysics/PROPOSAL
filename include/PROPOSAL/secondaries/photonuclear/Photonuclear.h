#pragma once
#include <stdexcept>

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/RegisteredInDefault.h"

namespace PROPOSAL {
namespace secondaries {
    struct Photonuclear : public Parametrization,
                          public RegisteredInDefault<Photonuclear> {
        int primary_particle_type;

        Photonuclear(ParticleDef p, Medium)
            : primary_particle_type(p.particle_type){};

        static constexpr size_t n_rnd = 0;
        size_t RequiredRandomNumbers() const noexcept { return n_rnd; }

        static constexpr InteractionType type
            = PROPOSAL::InteractionType::Photonuclear;
        InteractionType GetInteractionType() const noexcept { return type; }

        vector<Loss::secondary_t> CalculateSecondaries(double primary_energy,
            Loss::secondary_t loss, const Component&, vector<double>)
        {
            std::cout << "a hadronic interaction modell must be called."
                      << std::endl;
            std::get<Loss::TYPE>(loss) = primary_particle_type;
            std::get<Loss::ENERGY>(loss)
                = primary_energy - std::get<Loss::ENERGY>(loss);
            return std::vector<Loss::secondary_t>{loss};
        };
    };
} // namespace secondaries
} // namespace PROPOSAL
