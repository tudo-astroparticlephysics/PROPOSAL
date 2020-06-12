#pragma once

#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/particle/Particle.h"

using std::tuple;

namespace PROPOSAL {
namespace secondaries {
    struct Compton {
        Compton() = default;
        virtual ~Compton() = default;

        static constexpr InteractionType type
            = PROPOSAL::InteractionType::Compton;

        virtual double CalculateRho(double, double) = 0;
        virtual tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double)
            = 0;
        virtual tuple<double, double> CalculateEnergy(double, double) = 0;
    };
} // namespace secondaries
} // namespace PROPOSAL
