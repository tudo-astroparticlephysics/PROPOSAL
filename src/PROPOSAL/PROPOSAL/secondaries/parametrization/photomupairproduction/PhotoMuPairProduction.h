#pragma once

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/parametrization/Parametrization.h"

namespace PROPOSAL {
    namespace secondaries {
        struct PhotoMuPairProduction : public secondaries::Parametrization {
            PhotoMuPairProduction() = default;
            virtual ~PhotoMuPairProduction() = default;

            static constexpr InteractionType type
                = PROPOSAL::InteractionType::PhotoMuPair;
            InteractionType GetInteractionType() const noexcept { return type; };

            virtual double Calculatex(double, double, const Component&) = 0;
            virtual std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
                    const Vector3D&, double, const Component&, double) = 0;
            virtual std::tuple<double, double> CalculateEnergy(double, double)
            = 0;
        };
    } // namespace secondaries
} // namespace PROPOSAL
