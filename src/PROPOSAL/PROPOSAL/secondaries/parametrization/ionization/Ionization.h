#pragma once

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/parametrization/Parametrization.h"

namespace PROPOSAL {
namespace secondaries {
    struct Ionization : public secondaries::Parametrization {
        Ionization() = default;
        virtual ~Ionization() = default;

        static constexpr InteractionType type = PROPOSAL::InteractionType::Ioniz;
        InteractionType GetInteractionType() const noexcept { return type; };

        virtual std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
            const Vector3D&, double, double, double)
            = 0;
        virtual std::tuple<double, double> CalculateEnergy(double, double) = 0;
    };
} // namespace secondaries
} // namespace PROPOSAL
