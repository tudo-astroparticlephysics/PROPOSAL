
#pragma once

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/secondaries/ionization/Ioniaation.h"

/* #include <unordered_map> */

using PROPOSAL::Components::Component;
/* using std::unordered_map; */
using std::array;
using std::vector;

namespace PROPOSAL {
namespace secondaries {
    struct NaivIonization : public secondaries::Ionization {
        int primary_particle_type;

        static constexpr int n_rnd = 2;

        NaivIonization(const ParticleDef&);

        double CalculateRho(double, double) final;
        tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double) final;
        tuple<double, double> CalculateEnergy(double, double) final;

        virtual vector<Loss::secondary_t> CalculateSecondaries(
            Loss::secondary_t, array<double, n_rnd>);
    };
} // namespace secondaries
} // namespace PROPOSAL
