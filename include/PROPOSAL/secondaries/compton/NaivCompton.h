#pragma once

#include "PROPOSAL/secondaries/compton/Compton.h"
#include "PROPOSAL/secondaries/DefaultFactory.h"

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
        enum { GAMMA, ELECTRON };
        tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double) final;
        tuple<double, double> CalculateEnergy(double, double) final;

        size_t RequiredRandomNumbers() const noexcept { return n_rnd; }
        vector<Loss::secondary_t> CalculateSecondaries(
            double, Loss::secondary_t, const Component&, vector<double>);
    };
} // namespace secondaries
} // namespace PROPOSAL
