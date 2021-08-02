#pragma once

#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/NaivBremsstrahlung.h"

namespace PROPOSAL {
    namespace secondaries {
        class BremsKochMotz
                : public secondaries::NaivBremsstrahlung {
            static constexpr int n_rnd = 2;

            double CalculatePhiCandidate(double, double);
            double g(double, double, double, double);
            double N(double, double, double);
        public:
            BremsKochMotz() = delete;
            BremsKochMotz(ParticleDef p, Medium m)
                    : secondaries::NaivBremsstrahlung(p, m) {};

            size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }

            std::pair<Cartesian3D, Cartesian3D> CalculateDirections(
                    const Vector3D&, double, double, const Component&,
                    std::vector<double>&) final;
        };
    } // namespace secondaries
} // namespace PROPOSAL
