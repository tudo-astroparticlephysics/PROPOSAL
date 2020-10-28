#pragma once

#include "PROPOSAL/crosssection/parametrization/MupairProduction.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/secondaries/mupairproduction/MupairProduction.h"
#include "PROPOSAL/secondaries/RegisteredInDefault.h"
/* #include <unordered_map> */

using PROPOSAL::Components::Component;
/* using std::unordered_map; */
using std::array;
using std::unique_ptr;

namespace PROPOSAL {
namespace secondaries {
    class KelnerKokoulinPetrukhinMupairProduction
        : public secondaries::MupairProduction,
          public RegisteredInDefault<KelnerKokoulinPetrukhinMupairProduction> {

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
        tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double) final;
        tuple<double, double> CalculateEnergy(double, double) final;

        size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }
        vector<ParticleState> CalculateSecondaries(
                StochasticLoss, const Component&, vector<double>&) final;
    };
} // namespace secondaries
} // namespace PROPOSAL
