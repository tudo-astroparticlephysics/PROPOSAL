#pragma once

#include "PROPOSAL/secondaries/Parametrization.h"
#include "PROPOSAL/secondaries/DefaultFactory.h"

using std::tuple;
using std::vector;

namespace PROPOSAL {
namespace secondaries {
    struct MupairProduction : public secondaries::Parametrization  {
        MupairProduction() = default;
        virtual ~MupairProduction() = default;

        static InteractionType GetInteractionType() {
             return PROPOSAL::InteractionType::MuPair;
        };

        virtual double CalculateRho(double, double, const Component&, double) = 0;
        virtual tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double)
            = 0;
        virtual tuple<double, double> CalculateEnergy(double, double) = 0;
    };
} // namespace secondaries
} // namespace PROPOSAL
