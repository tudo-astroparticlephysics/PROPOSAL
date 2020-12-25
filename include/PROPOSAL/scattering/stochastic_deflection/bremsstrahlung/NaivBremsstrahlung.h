#pragma once

#include "PROPOSAL/scattering/stochastic_deflection/bremsstrahlung/Bremsstrahlung.h"

namespace PROPOSAL {
namespace stochastic_deflection {
    class NaivBremsstrahlung : public Bremsstrahlung,
                               public DefaultDeflection<NaivBremsstrahlung> {

        static constexpr int n_rnd = 1;

        std::unique_ptr<Parametrization> clone() const final
        {
            return std::unique_ptr<Parametrization>(
                std::make_unique<NaivBremsstrahlung>(*this));
        }

    public:
        NaivBremsstrahlung(ParticleDef, Medium) {};

        size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }

        std::array<double, 2> CalculateStochasticDeflection(
            double e_i, double e_f, std::vector<double> const& rnd) const final;
    };
} // namespace stochastic_deflection
} // namespace PROPOSAL
