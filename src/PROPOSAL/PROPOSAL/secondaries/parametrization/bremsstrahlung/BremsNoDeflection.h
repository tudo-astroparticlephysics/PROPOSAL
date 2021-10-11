#pragma once

#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/Bremsstrahlung.h"

namespace PROPOSAL {
namespace secondaries {
    class BremsNoDeflection
        : public secondaries::Bremsstrahlung {
        static constexpr int n_rnd = 0;
    public:
        BremsNoDeflection() = delete;
        BremsNoDeflection(const ParticleDef& p, const Medium&)
            : secondaries::Bremsstrahlung(p) {};

        size_t RequiredRandomNumbers() const noexcept override { return n_rnd; }

        std::pair<Cartesian3D, Cartesian3D> CalculateDirections(
                const Vector3D&, double, double,
                const Component&, std::vector<double>&) final;
    };
} // namespace secondaries
} // namespace PROPOSAL
