#pragma once

#include "PROPOSAL/secondaries/parametrization/compton/Compton.h"
#include "PROPOSAL/secondaries/RegisteredInDefault.h"

#include <vector>

using std::array;
using std::vector;

namespace PROPOSAL {
namespace secondaries {
    struct NaivCompton : public Compton,
                         public RegisteredInDefault<NaivCompton> {
        NaivCompton() = default;
        NaivCompton(ParticleDef, Medium) {};

        static constexpr int n_rnd = 1;

        double CalculateRho(double, double) final;
        enum { GAMMA, EMINUS };
        std::tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double) final;
        std::tuple<double, double> CalculateEnergy(double, double) final;

        size_t RequiredRandomNumbers() const noexcept { return n_rnd; }
        vector<ParticleState> CalculateSecondaries(StochasticLoss, const Component&,
                                                   vector<double>&);
    };
} // namespace secondaries
} // namespace PROPOSAL
