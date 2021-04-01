#pragma once

#include "PROPOSAL/secondaries/parametrization/Parametrization.h"

namespace PROPOSAL {
namespace secondaries {
    struct MupairProduction : public secondaries::Parametrization {
        MupairProduction() = default;
        virtual ~MupairProduction() = default;

        static constexpr InteractionType type = PROPOSAL::InteractionType::MuPair;
        InteractionType GetInteractionType() const noexcept { return type; };

        virtual double CalculateRho(double, double, const Component&, double,
                                    double) = 0;
        virtual std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
            const Vector3D&, double, double, double)
            = 0;
        virtual std::tuple<double, double> CalculateEnergy(double, double) = 0;
    };
} // namespace secondaries
} // namespace PROPOSAL
