#pragma once

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/Parametrization.h"

namespace PROPOSAL {
namespace secondaries {
    struct Ionization : public secondaries::Parametrization {
        Ionization() = default;
        virtual ~Ionization() = default;

        static InteractionType GetInteractionType() {
             return PROPOSAL::InteractionType::Ioniz;
        };

        virtual double CalculateRho(double, double) = 0;
        virtual tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double) = 0;
        virtual tuple<double, double> CalculateEnergy(double, double) = 0;
    };
} // namespace secondaries
} // namespace PROPOSAL
