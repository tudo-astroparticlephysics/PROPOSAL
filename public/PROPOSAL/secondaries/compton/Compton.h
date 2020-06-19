#pragma once

#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/Parametrization.h"

using std::tuple;

namespace PROPOSAL {
namespace secondaries {
    struct Compton : public secondaries::Parametrization {
        Compton() = default;
        virtual ~Compton() = default;

        static InteractionType GetInteractionType()
        {
            return PROPOSAL::InteractionType::Compton;
        };

        virtual double CalculateRho(double, double) = 0;
        virtual tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double)
            = 0;
        virtual tuple<double, double> CalculateEnergy(double, double) = 0;
    };
} // namespace secondaries
} // namespace PROPOSAL
