#include "PROPOSAL/scattering/stochastic_deflection/ionization/NaivIonization.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

DirectionChangeAngular
stochastic_deflection::NaivIonization::CalculateStochasticDeflection(
        double e_i, double e_f, std::vector<double> const& rnd) const
{
    auto p_i = std::sqrt((e_i + mass) * (e_i - mass));
    auto p_f = std::sqrt((e_f + mass) * (e_f - mass));
    auto costheta = ((e_i + ME) * e_f - e_i * ME - mass * mass) / (p_i * p_f);

    return DirectionChangeAngular { std::acos(costheta), 2 * PI * rnd[0] };
}