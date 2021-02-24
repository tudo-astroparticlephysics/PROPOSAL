#pragma once

#include "PROPOSAL/crosssection/parametrization/EpairProduction.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/secondaries/parametrization/epairproduction/EpairProduction.h"

namespace PROPOSAL {
    namespace secondaries {
        class KelnerKokoulinPetrukhinEpairProduction
                : public secondaries::EpairProduction,
                  public DefaultSecondaries<KelnerKokoulinPetrukhinEpairProduction> {

            crosssection::EpairKelnerKokoulinPetrukhin param;
            Integral integral;
            ParticleDef p_def;
            static constexpr int n_rnd = 3;

            double CalculateRho(double, double, const Component&, double, double);
            std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
                    const Vector3D&, double, double, double);
            std::tuple<double, double> CalculateEnergy(double, double);

        public:
            KelnerKokoulinPetrukhinEpairProduction(
              const ParticleDef& p, const Medium&) : param(false), p_def(p) {}
            //TODO: set lpm to true when possible

            size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }
            std::vector<ParticleState> CalculateSecondaries(StochasticLoss, const Component&,
                                                       std::vector<double>&) final;
        };
    } // namespace secondaries
} // namespace PROPOSAL
