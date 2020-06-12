#pragma once

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/medium/Components.h"

#include <vector>
#include <array>

using std::tuple;
using std::array;
using PROPOSAL::Components::Component;

namespace PROPOSAL {
namespace secondaries {
    struct PhotopairProduction {
        PhotopairProduction() = default;
        virtual ~PhotopairProduction() = default;

        static constexpr InteractionType type
            = PROPOSAL::InteractionType::Photopair;

        virtual double CalculateRho(double, double) = 0;
        virtual tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, const Component&, array<double,3>)
            = 0;
        virtual tuple<double, double> CalculateEnergy(double, double) = 0;
    };
} // namespace secondaries
} // namespace PROPOSAL
