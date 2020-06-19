
#pragma once

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/secondaries/DefaultFactory.h"

/* #include <unordered_map> */

using PROPOSAL::Components::Component;
/* using std::unordered_map; */
using std::array;
using std::vector;

namespace PROPOSAL {
namespace secondaries {
    struct NaivEpairProduction
        : public secondaries::EpairProduction,
          public RegisteredInDefault<NaivEpairProduction> {
        // unordered_map<Component*, unique_ptr<Interpolant>> dndx;

        static constexpr int n_rnd = 2;

        NaivEpairProduction(ParticleDef, Medium);

        double CalculateRho(double, double) final;
        tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double) final;
        tuple<double, double> CalculateEnergy(double, double) final;

        virtual vector<Loss::secondary_t> CalculateSecondaries(
            Loss::secondary_t, array<double, n_rnd>);
    };
} // namespace secondaries
} // namespace PROPOSAL
