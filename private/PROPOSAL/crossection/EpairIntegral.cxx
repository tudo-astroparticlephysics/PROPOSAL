
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/EpairIntegral.h"
#include "PROPOSAL/crossection/parametrization/EpairProduction.h"
#include "PROPOSAL/medium/Medium.h"

using namespace PROPOSAL;

double EpairIntegral::FunctionToDEdxIntegralReverse(double energy, double v)
{
    return (1 - v) * parametrization_->DifferentialCrossSection(energy, v);
}

double EpairIntegral::dedx_integral(double energy)
{
    auto v_cut = GetEnergyCut(energy);
    auto v_min = parametrization_->GetKinematicLimits(energy).vMin;

    double r1 = 0.8;
    double rUp = v_cut * (1 - HALF_PRECISION);
    bool rflag = false;

    if (r1 < rUp) {
        if (2 * parametrization_->FunctionToDEdxIntegral(energy, r1)
            < parametrization_->FunctionToDEdxIntegral(energy, rUp)) {
            rflag = true;
        }
    }

    auto func = [&, energy](double v) {
        return parametrization_->FunctionToDEdxIntegral(energy, v);
    };

    if (rflag) {
        if (r1 > v_cut) {
            r1 = v_cut;
        }

        if (r1 < v_min) {
            r1 = v_min;
        }

        auto sum = integral_.Integrate(v_min, r1, func, 4);
        double r2 = std::max(1 - v_cut, COMPUTER_PRECISION);

        if (r2 > 1 - r1) {
            r2 = 1 - r1;
        }

        auto func_reverse = [&, energy](double v) {
            return FunctionToDEdxIntegralReverse(energy, v);
        };
        sum += integral_.Integrate(1 - v_cut, r2, func_reverse, 2)
            + integral_.Integrate(r2, 1 - r1, func_reverse, 4);
        return energy * sum;
    }
    return energy * integral_.Integrate(v_min, v_cut, func, 4);
}
