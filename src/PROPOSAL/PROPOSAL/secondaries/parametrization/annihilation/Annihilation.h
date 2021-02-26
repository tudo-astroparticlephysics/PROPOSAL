#pragma once

#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/parametrization/Parametrization.h"

namespace PROPOSAL {
namespace secondaries {
    struct Annihilation : public secondaries::Parametrization {
        Annihilation() = default;
        virtual ~Annihilation() = default;

        virtual double CalculateRho(double, double, const Component&) = 0;
        virtual std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
            const Vector3D&, double, double, double)
            = 0;
        virtual std::tuple<double, double> CalculateEnergy(double, double) = 0;

        static constexpr InteractionType type
            = PROPOSAL::InteractionType::Annihilation;

        InteractionType GetInteractionType() const noexcept final { return type; };
    };
} // namespace secondaries
} // namespace PROPOSAL
