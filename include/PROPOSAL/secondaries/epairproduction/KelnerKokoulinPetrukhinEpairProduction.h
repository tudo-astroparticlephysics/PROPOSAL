#pragma once

#include "PROPOSAL/crosssection/parametrization/EpairProduction.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/secondaries/epairproduction/EpairProduction.h"
#include "PROPOSAL/secondaries/DefaultFactory.h"

using PROPOSAL::Components::Component;
using std::array;
using std::unique_ptr;

namespace PROPOSAL {
    namespace secondaries {
        class KelnerKokoulinPetrukhinEpairProduction
                : public secondaries::EpairProduction,
                  public RegisteredInDefault<KelnerKokoulinPetrukhinEpairProduction> {

            crosssection::EpairKelnerKokoulinPetrukhin param;
            Integral integral;
            ParticleDef p_def;
            static constexpr int n_rnd = 2;

            double CalculateRho(double, double, const Component&, double);
            tuple<Vector3D, Vector3D> CalculateDirections(
                    Vector3D, double, double, double);
            tuple<double, double> CalculateEnergy(double, double);

        public:
            KelnerKokoulinPetrukhinEpairProduction(
              const ParticleDef& p, const Medium&) : param(false), p_def(p) {}
            //TODO: set lpm to true when possible

            size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }
            vector<Loss::secondary_t> CalculateSecondaries(
                    double, Loss::secondary_t, const Component&, vector<double>) final;
        };
    } // namespace secondaries
} // namespace PROPOSAL
