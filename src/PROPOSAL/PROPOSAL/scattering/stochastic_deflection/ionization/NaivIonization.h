#pragma once

#include "PROPOSAL/scattering/stochastic_deflection/ionization/Ionization.h"

namespace PROPOSAL {
    namespace stochastic_deflection {
        class NaivIonization : public Ionization,
                               public DefaultDeflection<NaivIonization> {

            static constexpr int n_rnd = 1;
            double mass;

            std::unique_ptr<Parametrization> clone() const final
            {
                return std::unique_ptr<Parametrization>(
                        std::make_unique<NaivIonization>(*this));
            }

        public:
            NaivIonization(ParticleDef p_def, Medium) : mass(p_def.mass) {};

            size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }

            std::array<double, 2> CalculateStochasticDeflection(
                    double e_i, double e_f, std::vector<double> const& rnd) const final;
        };
    } // namespace stochastic_deflection
} // namespace PROPOSAL
