#pragma once

#include "PROPOSAL/crosssection/parametrization/EpairProduction.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/secondaries/parametrization/epairproduction/EpairProduction.h"
#include "PROPOSAL/secondaries/RegisteredInDefault.h"

using PROPOSAL::Components::Component;
using std::array;
using std::unique_ptr;
using std::tuple;

namespace PROPOSAL {
    namespace secondaries {
        class KelnerKokoulinPetrukhinEpairProduction
                : public secondaries::EpairProduction,
                  public RegisteredInDefault<KelnerKokoulinPetrukhinEpairProduction> {

            crosssection::EpairKelnerKokoulinPetrukhin param;
            Integral integral;
            ParticleDef p_def;
            static constexpr int n_rnd = 3;

            double CalculateRho(double, double, const Component&, double, double);
            tuple<Vector3D, Vector3D> CalculateDirections(
                    Vector3D, double, double, double);
            tuple<double, double> CalculateEnergy(double, double);

        public:
            KelnerKokoulinPetrukhinEpairProduction(
              const ParticleDef& p, const Medium&) : param(false), p_def(p) {}
            //TODO: set lpm to true when possible

            size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }
            vector<ParticleState> CalculateSecondaries(StochasticLoss, const Component&,
                                                       vector<double>&) final;
        };
    } // namespace secondaries
} // namespace PROPOSAL
