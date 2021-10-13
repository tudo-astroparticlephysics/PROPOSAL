#pragma once

#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/Bremsstrahlung.h"

namespace PROPOSAL {
    namespace secondaries {
        class BremsKochMotz
                : public secondaries::Bremsstrahlung {
            static constexpr int n_rnd = 2;
            double primary_lepton_mass;

            double CalculatePhiCandidate(double, double);
            double g(double, double, double, double);
            double N(double, double, double);
        public:
            BremsKochMotz() = delete;
            BremsKochMotz(const ParticleDef& p, const Medium& m)
                    : secondaries::Bremsstrahlung(p),
                      primary_lepton_mass(p.mass) {};

            size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }

            std::pair<Cartesian3D, Cartesian3D> CalculateDirections(
                    const Vector3D&, double, double, const Component&,
                    std::vector<double>&) final;
        };
    } // namespace secondaries
} // namespace PROPOSAL
