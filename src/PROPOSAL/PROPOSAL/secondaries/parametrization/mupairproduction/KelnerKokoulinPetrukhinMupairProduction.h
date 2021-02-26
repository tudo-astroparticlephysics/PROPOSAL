#pragma once

#include "PROPOSAL/crosssection/parametrization/MupairProduction.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/secondaries/parametrization/mupairproduction/MupairProduction.h"

namespace PROPOSAL {
namespace secondaries {
    class KelnerKokoulinPetrukhinMupairProduction
        : public secondaries::MupairProduction,
          public DefaultSecondaries<KelnerKokoulinPetrukhinMupairProduction> {

        crosssection::MupairKelnerKokoulinPetrukhin param;
        Integral integral;
        ParticleDef p_def;

        static constexpr int n_rnd = 3;

    public:
        KelnerKokoulinPetrukhinMupairProduction(
            const ParticleDef& p, const Medium&)
            : p_def(p)
        {
        }

        double CalculateRho(double, double, const Component&, double, double) final;
        std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
            const Vector3D&, double, double, double) final;
        std::tuple<double, double> CalculateEnergy(double, double) final;

        size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }
        std::vector<ParticleState> CalculateSecondaries(
                StochasticLoss, const Component&, std::vector<double>&) final;
    };
} // namespace secondaries
} // namespace PROPOSAL
