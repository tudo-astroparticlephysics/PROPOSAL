#pragma once

#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/Bremsstrahlung.h"

namespace PROPOSAL {
namespace secondaries {
    class NaivBremsstrahlung
        : public secondaries::Bremsstrahlung {
        static constexpr int n_rnd = 0;
        const int primary_lepton_type;
    protected:
        double primary_lepton_mass;
    public:
        NaivBremsstrahlung() = delete;
        NaivBremsstrahlung(ParticleDef p, Medium)
            : primary_lepton_type(p.particle_type),
              primary_lepton_mass(p.mass) {};

        size_t RequiredRandomNumbers() const noexcept override { return n_rnd; }
        std::vector<ParticleState> CalculateSecondaries(
            StochasticLoss, const Component&, std::vector<double>&) override;

        virtual std::pair<Cartesian3D, Cartesian3D> CalculateDirections(
                const Vector3D&, double, double,
                const Component&, std::vector<double>&);
    };
} // namespace secondaries
} // namespace PROPOSAL
