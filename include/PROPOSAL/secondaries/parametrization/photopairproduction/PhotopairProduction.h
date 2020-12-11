#pragma once

#include "PROPOSAL/medium/Components.h"
/* #include "PROPOSAL/particle/Particle.h" */
#include "PROPOSAL/secondaries/parametrization/Parametrization.h"

#include <array>
#include <vector>

using PROPOSAL::Components::Component;
using std::array;
using std::tuple;

namespace PROPOSAL {
namespace secondaries {
    struct PhotopairProduction : public secondaries::Parametrization {
        PhotopairProduction() = default;
        virtual ~PhotopairProduction() = default;

        static constexpr InteractionType type
            = PROPOSAL::InteractionType::Photopair;
        InteractionType GetInteractionType() const noexcept { return type; };

        virtual double CalculateRho(double, double, const Component&) = 0;
        virtual tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, const Component&, vector<double>)
            = 0;
        virtual tuple<double, double> CalculateEnergy(double, double, double)
            = 0;
    };
} // namespace secondaries
} // namespace PROPOSAL
