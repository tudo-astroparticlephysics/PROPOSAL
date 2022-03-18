#pragma once

#include "PROPOSAL/scattering/stochastic_deflection/epairproduction/EpairProduction.h"

namespace PROPOSAL {
    namespace stochastic_deflection {
        class EpairGinneken : public EpairProduction,
                public DefaultDeflection<EpairGinneken> {

            static constexpr int n_rnd = 2;
            double mass;

            std::unique_ptr<Parametrization> clone() const final
            {
                return std::make_unique<EpairGinneken>(*this);
            }

        public:
            EpairGinneken(const ParticleDef& p_def, const Medium&)
                : mass(p_def.mass){};

            size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }

            UnitSphericalVector CalculateStochasticDeflection(
                    double e_i, double e_f, std::vector<double> const& rnd, size_t) const final;
        };
    } // namespace stochastic_deflection
} // namespace PROPOSAL

// S.R. Kel'ner, Sov. J. Nucl. Phys.5 (1967)778.