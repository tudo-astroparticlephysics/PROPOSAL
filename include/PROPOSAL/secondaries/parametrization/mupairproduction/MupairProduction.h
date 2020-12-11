#pragma once

#include "PROPOSAL/secondaries/parametrization/Parametrization.h"

using std::tuple;
using std::vector;

namespace PROPOSAL {
namespace secondaries {
    struct MupairProduction : public secondaries::Parametrization {
        MupairProduction() = default;
        virtual ~MupairProduction() = default;

        static constexpr InteractionType type = PROPOSAL::InteractionType::MuPair;
        InteractionType GetInteractionType() const noexcept { return type; };

        virtual double CalculateRho(double, double, const Component&, double,
                                    double) = 0;
        virtual tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double)
            = 0;
        virtual tuple<double, double> CalculateEnergy(double, double) = 0;
    };
} // namespace secondaries
} // namespace PROPOSAL
