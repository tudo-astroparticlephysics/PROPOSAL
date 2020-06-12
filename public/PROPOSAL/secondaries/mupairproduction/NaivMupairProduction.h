
#pragma once

#include "PROPOSAL/secondaries/mupairproduction/MupairProduction.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/math/Integral.h"
/* #include <unordered_map> */

using PROPOSAL::Components::Component;
/* using std::unordered_map; */
using std::array;

namespace PROPOSAL {
namespace secondaries {
    class NaivMupairProduction : public secondaries::MupairProduction {
        Integral integral;
        double(*rho_integrand)(double, double, double);

        static constexpr int n_rnd = 2;

        NaivMupairProduction() = default;

        double CalculateRho(double, double, double) final;
        tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, double) final;
        tuple<double, double> CalculateEnergy(double, double) final;

        virtual vector<Loss::secondary_t> CalculateSecondaries(
            Loss::secondary_t, array<double, n_rnd>);
    };
} // namespace secondaries
} // namespace PROPOSAL
