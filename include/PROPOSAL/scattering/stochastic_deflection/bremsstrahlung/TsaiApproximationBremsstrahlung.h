#pragma once

#include "PROPOSAL/scattering/stochastic_deflection/bremsstrahlung/Bremsstrahlung.h"

namespace PROPOSAL {
namespace stochastic_deflection {
    class TsaiApproximationBremsstrahlung : public Bremsstrahlung,
                               public DefaultDeflection<TsaiApproximationBremsstrahlung> {

        static constexpr int n_rnd = 2;
        double mass;

        std::unique_ptr<Parametrization> clone() const final
        {
            return std::unique_ptr<Parametrization>(
                std::make_unique<TsaiApproximationBremsstrahlung>(*this));
        }

    public:
        TsaiApproximationBremsstrahlung(const ParticleDef& p_def, const Medium&)
            : mass(p_def.mass) {};

        size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }

        DirectionChangeAngular CalculateStochasticDeflection(
            double e_i, double e_f, std::vector<double> const& rnd) const final;
    };
} // namespace stochastic_deflection
} // namespace PROPOSAL

// SR Kelner, RP Kokoulin, and AA Petrukhin. About cross section for high-energy muon bremsstrahlung. Technical Report, MEphI, 1995. Preprint MEPhI 024-95, Moscow, 1995, CERN SCAN-9510048.
// R.P. Kokoulin S.R. Kelner and A.A. Petrukhin. Bremsstrahlung from muons scattered by atomic elec- trons. Physics of Atomic Nuclei, 60:576–583, April 1997.
