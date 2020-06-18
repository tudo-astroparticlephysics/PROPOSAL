#pragma once

#include "PROPOSAL/crossection/parametrization/MupairProduction.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/secondaries/mupairproduction/MupairProduction.h"
/* #include <unordered_map> */

using PROPOSAL::Components::Component;
/* using std::unordered_map; */
using std::array;
using std::unique_ptr;

namespace PROPOSAL {
namespace secondaries {
    class KelnerKokoulinPetrukhinMupairProduction
        : public secondaries::MupairProduction {
        crosssection::MupairKelnerKokoulinPetrukhin param;
        Integral integral;
        const ParticleDef p_def;

        static constexpr int n_rnd = 2;

    public:
        KelnerKokoulinPetrukhinMupairProduction(ParticleDef);

        double CalculateRho(double, double, const Component&, double) final;
        tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double) final;
        tuple<double, double> CalculateEnergy(double, double) final;

        size_t RequiredRandomNumbers() final { return n_rnd; }
        virtual vector<Loss::secondary_t> CalculateSecondaries(
            double, Loss::secondary_t, const Component&, vector<double>);
    };
} // namespace secondaries
} // namespace PROPOSAL
