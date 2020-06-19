#pragma once

#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/Parametrization.h"

using std::tuple;

namespace PROPOSAL {
namespace secondaries {
    struct Annihilation : public secondaries::Parametrization {
        Annihilation() = default;
        ~Annihilation() = default;

        virtual double CalculateRho(double, double, const Component&) = 0;
        virtual tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double)
            = 0;
        virtual tuple<double, double> CalculateEnergy(double, double) = 0;

        static InteractionType GetInteractionType()
        {
            return PROPOSAL::InteractionType::Annihilation;
        };
    };
} // namespace secondaries
} // namespace PROPOSAL
