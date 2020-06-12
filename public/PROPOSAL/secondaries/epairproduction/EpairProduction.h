#pragma once

#include "PROPOSAL/particle/Particle.h"

namespace PROPOSAL {
namespace secondaries {
    struct EpairProduction {
        EpairProduction() = default;
        virtual ~EpairProduction() = default;

        static constexpr InteractionType type
            = PROPOSAL::InteractionType::Epair;

        virtual double CalculateRho(double, double) = 0;
        virtual tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double)
            = 0;
        virtual tuple<double, double> CalculateEnergy(double, double) = 0;
    };
} // namespace secondaries
} // namespace PROPOSAL
