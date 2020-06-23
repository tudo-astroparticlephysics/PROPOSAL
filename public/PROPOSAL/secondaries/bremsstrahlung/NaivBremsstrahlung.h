#pragma once

#include "PROPOSAL/secondaries/bremsstrahlung/Bremsstrahlung.h"
#include "PROPOSAL/secondaries/DefaultFactory.h"

using std::array;
using std::vector;

namespace PROPOSAL {
namespace secondaries {
    class NaivBremsstrahlung : public secondaries::Bremsstrahlung,
                               public RegisteredInDefault<NaivBremsstrahlung> {
        static constexpr int n_rnd = 0;

    public:
        NaivBremsstrahlung() = default;
        NaivBremsstrahlung(ParticleDef, Medium) {};

        size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }
        vector<Loss::secondary_t> CalculateSecondaries(
            double, Loss::secondary_t, const Component&, vector<double>);
    };
} // namespace secondaries
} // namespace PROPOSAL
