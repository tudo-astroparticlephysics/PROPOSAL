#pragma once

#include "PROPOSAL/scattering/stochastic_deflection/pairProd/PairProduction.h"

namespace PROPOSAL {
    namespace stochastic_deflection {
        class KelnerPairProduction : public PairProduction,
                               public DefaultDeflection<KelnerPairProduction> {

            static constexpr int n_rnd = 2;

            std::unique_ptr<Parametrization> clone() const final
            {
                return std::unique_ptr<Parametrization>(
                        std::make_unique<KelnerPairProduction>(*this));
            }

        public:
            KelnerPairProduction(ParticleDef, Medium) {};

            size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }

            std::array<double, 2> CalculateStochasticDeflection(
                    double e_i, double e_f, std::vector<double> const& rnd) const final;
        };
    } // namespace stochastic_deflection
} // namespace PROPOSAL

// S.R. Kel'ner, Sov. J. Nucl. Phys.5 (1967)778.