#pragma once

#include "PROPOSAL/secondaries/compton/Compton.h"

#include <vector>

using std::vector;
using std::array;

namespace PROPOSAL {
namespace secondaries {
    struct NaivCompton : public Compton{
        NaivCompton() = default;

        static constexpr int n_rnd = 1;

        double CalculateRho(double, double) final;
        enum {GAMMA, ELECTRON};
        tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double) final;
        tuple<double, double> CalculateEnergy(double, double) final;

        vector<Loss::secondary_t> CalculateSecondaries(
            Loss::secondary_t, array<double, n_rnd>);
    };
} // namespace secondaries
} // namespace PROPOSAL
