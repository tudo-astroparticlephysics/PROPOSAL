#pragma once

#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/parametrization/Parametrization.h"

using std::tuple;

namespace PROPOSAL {
namespace secondaries {
    struct Annihilation : public secondaries::Parametrization {
        Annihilation() = default;
        virtual ~Annihilation() = default;

        virtual double CalculateRho(double, double, const Component&) = 0;
        virtual tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double)
            = 0;
        virtual tuple<double, double> CalculateEnergy(double, double) = 0;

        static constexpr InteractionType type
            = PROPOSAL::InteractionType::Annihilation;

        InteractionType GetInteractionType() const noexcept final { return type; };
    };
} // namespace secondaries
} // namespace PROPOSAL
