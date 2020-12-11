#pragma once

#include "PROPOSAL/scattering/stochastic_deflection/RegisteredInDefault.h"
#include "PROPOSAL/scattering/stochastic_deflection/parametrization/bremsstrahlung/Bremsstrahlung.h"

namespace PROPOSAL {
namespace stochastic_deflection {
    class NaivBremsstrahlung : public stochastic_deflection::Bremsstrahlung,
                               public RegisteredInDefault<NaivBremsstrahlung> {

        static constexpr int n_rnd = 0;

    public:
        NaivBremsstrahlung() = delete;
        NaivBremsstrahlung(ParticleDef p, Medium);

        size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }

        std::array<double, 2> CalculateDeflection(
            StochasticLoss, const Component&, std::vector<double>&);
    };
} // namespace stochastic_deflection
} // namespace PROPOSAL
