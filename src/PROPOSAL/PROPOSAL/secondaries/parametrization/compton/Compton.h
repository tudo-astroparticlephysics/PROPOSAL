#pragma once

#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/parametrization/Parametrization.h"

namespace PROPOSAL {
namespace secondaries {
    struct Compton : public secondaries::Parametrization {
        Compton() = default;
        virtual ~Compton() = default;

        static constexpr InteractionType type = PROPOSAL::InteractionType::Compton;

        InteractionType GetInteractionType() const noexcept { return type; };

        virtual double CalculateRho(double, double) = 0;
        virtual std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
            const Vector3D&, double, double, double)
            = 0;
        virtual std::tuple<double, double> CalculateEnergy(double, double) = 0;
    };
} // namespace secondaries
} // namespace PROPOSAL
