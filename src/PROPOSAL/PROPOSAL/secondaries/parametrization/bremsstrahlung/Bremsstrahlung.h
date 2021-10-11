#pragma once

#include "PROPOSAL/secondaries/parametrization/Parametrization.h"

namespace PROPOSAL {
namespace secondaries {
    struct Bremsstrahlung : public Parametrization {
        const int primary_lepton_type;
    public:
        Bremsstrahlung() = delete;
        Bremsstrahlung(const ParticleDef& p)
            : primary_lepton_type(p.particle_type) {};
        virtual ~Bremsstrahlung() = default;

        static constexpr InteractionType type = InteractionType::Brems;

        std::vector<ParticleState> CalculateSecondaries(
                StochasticLoss, const Component&, std::vector<double>&) override;
        virtual std::pair<Cartesian3D, Cartesian3D> CalculateDirections(
                const Vector3D&, double, double,
                const Component&, std::vector<double>&) = 0;
        InteractionType GetInteractionType() const noexcept final { return type; };
    };
} // namespace secondaries
} // namespace PROPOSAL
