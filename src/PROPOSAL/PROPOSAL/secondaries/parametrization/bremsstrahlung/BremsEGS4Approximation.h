#pragma once

#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/NaivBremsstrahlung.h"

namespace PROPOSAL {
namespace secondaries {
    class BremsEGS4Approximation
            : public secondaries::NaivBremsstrahlung,
              public DefaultSecondaries<BremsEGS4Approximation> {
        static constexpr int n_rnd = 1;
        double primary_lepton_mass;

    public:
        BremsEGS4Approximation() = delete;
        BremsEGS4Approximation(ParticleDef p, Medium m)
            : secondaries::NaivBremsstrahlung(p, m),
            primary_lepton_mass(p.mass) {};

        size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }

        std::pair<Cartesian3D, Cartesian3D> CalculateDirections(
                const Vector3D&, double, double, const Component&,
                std::vector<double>&) final;
    };
} // namespace secondaries
} // namespace PROPOSAL
