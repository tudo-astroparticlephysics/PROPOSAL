
#pragma once

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/secondaries/RegisteredInDefault.h"
#include "PROPOSAL/secondaries/ionization/Ionization.h"

/* #include <unordered_map> */

using PROPOSAL::Components::Component;
/* using std::unordered_map; */
using std::array;
using std::vector;

namespace PROPOSAL {
namespace secondaries {
    struct NaivIonization : public secondaries::Ionization,
                            public RegisteredInDefault<NaivIonization> {
        int primary_particle_type;

        static constexpr int n_rnd = 2;

        NaivIonization(const ParticleDef& p, const Medium&) :
            primary_particle_type(p.particle_type) {}

        tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double) final;
        tuple<double, double> CalculateEnergy(double, double) final;

        size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }
        vector<DynamicData> CalculateSecondaries(StochasticLoss, const Component&,
                                                 vector<double>&) final;
    };
} // namespace secondaries
} // namespace PROPOSAL
