#pragma once

#include "PROPOSAL/scattering/stochastic_deflection/bremsstrahlung/Bremsstrahlung.h"

namespace PROPOSAL {
namespace stochastic_deflection {
    class NaivBremsstrahlung : public Bremsstrahlung,
                               public DefaultDeflection<NaivBremsstrahlung> {

        static constexpr int n_rnd = 0;

    public:
        NaivBremsstrahlung(ParticleDef p, Medium) {};

        size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }

        std::array<double, 2> CalculateStochasticDeflection(
            double initial_energy, double final_energy, std::vector<double> const&) const final;
    };
} // namespace stochastic_deflection
} // namespace PROPOSAL
