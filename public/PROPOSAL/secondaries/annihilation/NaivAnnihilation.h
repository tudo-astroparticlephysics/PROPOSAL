#pragma once

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crossection/parametrization/Annihilation.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/secondaries/annihilation/Annihilation.h"

/* #include <unordered_map> */

using PROPOSAL::Components::Component;
/* using std::unordered_map; */
using std::vector;
using std::array;

namespace PROPOSAL {
namespace secondaries {
    struct NaivAnnihilation : public secondaries::Annihilation {
        // unordered_map<Component*, unique_ptr<Interpolant>> dndx;

        static constexpr int n_rnd = 2;

        NaivAnnihilation(const crosssection::Annihilation&, const ParticleDef&,
            const Medium&, const EnergyCutSettings&);

        double CalculateRho(double, double) final;
        tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double) final;
        tuple<double, double> CalculateEnergy(double, double) final;

        virtual vector<Loss::secondary_t> CalculateSecondaries(
            Loss::secondary_t, array<double, n_rnd>);
    };
} // namespace secondaries
} // namespace PROPOSAL
