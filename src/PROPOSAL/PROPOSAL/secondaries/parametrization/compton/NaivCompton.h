#pragma once

#include "PROPOSAL/secondaries/parametrization/compton/Compton.h"

#include <vector>

namespace PROPOSAL {
namespace secondaries {
    struct NaivCompton : public Compton,
                         public DefaultSecondaries<NaivCompton> {
        NaivCompton() = default;
        NaivCompton(ParticleDef, Medium) {};

        static constexpr int n_rnd = 1;

        double CalculateRho(double, double) final;
        enum { GAMMA, EMINUS };
        std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
            const Vector3D&, double, double, double) final;
        std::tuple<double, double> CalculateEnergy(double, double) final;

        size_t RequiredRandomNumbers() const noexcept { return n_rnd; }
        std::vector<ParticleState> CalculateSecondaries(StochasticLoss, const Component&,
                                                   std::vector<double>&);
    };
} // namespace secondaries
} // namespace PROPOSAL
