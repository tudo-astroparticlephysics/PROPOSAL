#pragma once

#include "PROPOSAL/secondaries/bremsstrahlung/Bremsstrahlung.h"

using std::array;
using std::vector;

namespace PROPOSAL {
namespace secondaries {
    struct NaivBremsstrahlung : public secondaries::Bremsstrahlung {
        static constexpr int n_rnd = 0;

        NaivBremsstrahlung() = default;

        vector<Loss::secondary_t> CalculateSecondaries(
            Loss::secondary_t, array<double, n_rnd>);
    };
} // namespace secondaries
} // namespace PROPOSAL
