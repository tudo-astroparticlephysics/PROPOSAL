#pragma once

#include "PROPOSAL/particle/Particle.h"

namespace PROPOSAL {
namespace secondaries {
    struct Ionization {
        Ionization() = default;
        virtual ~Ionization() = default;

        static constexpr InteractionType type
            = PROPOSAL::InteractionType::Ioniz;

        double CalculateRho(double, double) final;
        tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double) final;
        tuple<double, double> CalculateEnergy(double, double) final;
    };
} // namespace secondaries
} // namespace PROPOSAL
