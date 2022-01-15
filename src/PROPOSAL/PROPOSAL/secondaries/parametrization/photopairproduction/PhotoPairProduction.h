#pragma once

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/parametrization/Parametrization.h"

namespace PROPOSAL {
namespace secondaries {
    struct PhotoPairProduction : public secondaries::Parametrization {
        PhotoPairProduction() = default;
        virtual ~PhotoPairProduction() = default;

        static constexpr InteractionType type
            = PROPOSAL::InteractionType::Photopair;
        InteractionType GetInteractionType() const noexcept { return type; };

        virtual double CalculateRho(double, double, const Component&) = 0;
        virtual std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
            const Vector3D&, double, double, const Component&, double, double,
            double) = 0;
        virtual std::tuple<double, double> CalculateEnergy(double, double, double)
            = 0;
    };
} // namespace secondaries
} // namespace PROPOSAL
