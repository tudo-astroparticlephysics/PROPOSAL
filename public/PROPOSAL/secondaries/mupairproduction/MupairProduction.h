#pragma once

#include "PROPOSAL/particle/Particle.h"

using std::vector;
using std::tuple;

namespace PROPOSAL {
namespace secondaries {
    struct MupairProduction {
        MupairProduction() = default;
        virtual ~MupairProduction() = default;

        static constexpr InteractionType type
            = PROPOSAL::InteractionType::MuPair;

        virtual double CalculateRho(double, double, double) = 0;
        virtual tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double)
            = 0;
        virtual tuple<double, double> CalculateEnergy(double, double) = 0;
    };
} // namespace secondaries
} // namespace PROPOSAL
