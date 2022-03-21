#include "PROPOSAL/scattering/stochastic_deflection/ionization/IonizNaive.h"
#include "PROPOSAL/Constants.h"
#include <cmath>
using namespace PROPOSAL;

UnitSphericalVector
stochastic_deflection::IonizNaive::CalculateStochasticDeflection(
        double e_i, double e_f, std::vector<double> const& rnd, size_t) const
{
    auto p_i = std::sqrt((e_i + mass) * (e_i - mass));
    auto p_f = std::sqrt((e_f + mass) * (e_f - mass));
    auto costheta = ((e_i + ME) * e_f - e_i * ME - mass * mass) / (p_i * p_f);

    return UnitSphericalVector { std::acos(costheta), 2 * PI * rnd[0] };
}