#pragma once

#include "PROPOSAL/scattering/stochastic_deflection/nuclearInteraction/NuclearInteraction.h"


namespace PROPOSAL {
    namespace stochastic_deflection {
        class BorogPetrukhinNuclearInteraction : public NuclearInteraction,
                               public DefaultDeflection<BorogPetrukhinNuclearInteraction> {

            static constexpr int n_rnd = 2;
            double mass;

            std::unique_ptr<Parametrization> clone() const final
            {
                return std::unique_ptr<Parametrization>(
                        std::make_unique<BorogPetrukhinNuclearInteraction>(*this));
            }

        public:
            BorogPetrukhinNuclearInteraction(const ParticleDef& p_def, const Medium&)
                : mass(p_def.mass) {};

            size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }

            DirectionChangeAngular CalculateStochasticDeflection(
                    double e_i, double e_f, std::vector<double> const& rnd) const final;
        };
    } // namespace stochastic_deflection
} // namespace PROPOSAL

// V. V Borog, V. G. Kirillov-Ugryumov, and A. A Petrukhin. Inelastic Interactions of Muons with Fe Nuclei at Small q**2 in the Energy Region Above 200-GeV. Yad. Fiz., 25:85–93, 1977. Sov. J. Nucl. Phys.25, 46-51 (1977) - english translation.
// V.V. Borog and A.A. Petrukhin. The cross-section of the nuclear interaction of high energy muons. International Cosmic Ray Conference, 6:1949–1954, August 1975.