#pragma once

#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/Bremsstrahlung.h"

namespace PROPOSAL {
namespace secondaries {
    class BremsEGS4Approximation
            : public secondaries::Bremsstrahlung,
              public DefaultSecondaries<BremsEGS4Approximation> {
        static constexpr int n_rnd = 1;
        double primary_lepton_mass;

    public:
        BremsEGS4Approximation() = delete;
        BremsEGS4Approximation(const ParticleDef& p, const Medium&)
            : secondaries::Bremsstrahlung(p),
              primary_lepton_mass(p.mass) {};

        size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }

        std::pair<Cartesian3D, Cartesian3D> CalculateDirections(
                const Vector3D&, double, double, const Component&,
                std::vector<double>&) final;
    };
} // namespace secondaries
} // namespace PROPOSAL
