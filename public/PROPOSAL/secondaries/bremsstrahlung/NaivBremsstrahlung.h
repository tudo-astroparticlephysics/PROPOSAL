#pragma once

#include "PROPOSAL/secondaries/bremsstrahlung/Bremsstrahlung.h"

using std::array;
using std::vector;

namespace PROPOSAL {
namespace secondaries {
    class NaivBremsstrahlung : public secondaries::Bremsstrahlung {
        static constexpr int n_rnd = 0;

    public:
        NaivBremsstrahlung() = default;

        size_t RequiredRandomNumbers() { return n_rnd; }
        vector<Loss::secondary_t> CalculateSecondaries(
            double, Loss::secondary_t, const Component&, vector<double>);
    };
} // namespace secondaries
} // namespace PROPOSAL
