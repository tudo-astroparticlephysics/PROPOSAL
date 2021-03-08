#include "PROPOSAL/scattering/stochastic_deflection/nuclearInteraction/BorogPetrukhinNuclearInteraction.h"
#include "PROPOSAL/Constants.h"
#include <cmath>

using namespace PROPOSAL;

DirectionChangeAngular
stochastic_deflection::BorogPetrukhinNuclearInteraction::CalculateStochasticDeflection(
        double e_i, double e_f, std::vector<double> const& rnd) const
{
    auto m_0 = std::sqrt(0.4) * 1e3;
    auto epsilon = e_i - e_f; 
    auto y = epsilon / e_i; 
    auto t_max = 2 * MP * epsilon;
    auto t_min = mass * mass * y * y / (1.0 - y);
    auto t_1 = std::min(epsilon * epsilon, m_0 * m_0);
    auto t_p = (t_max * t_1) / ((t_max + t_1) * pow(((t_max * (t_min + t_1)) / (t_min * (t_max + t_1))), rnd[0]) - t_max); 
    auto sin2 = (t_p - t_min) / (4 * (e_i * e_f - mass * mass) - 2 * t_min);
    auto theta_muon = 2 * std::asin(std::sqrt(sin2));

    return DirectionChangeAngular { theta_muon, 2 * PI * rnd[1] };
}