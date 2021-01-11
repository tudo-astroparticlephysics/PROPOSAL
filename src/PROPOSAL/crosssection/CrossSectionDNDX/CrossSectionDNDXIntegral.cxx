#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"

using namespace PROPOSAL;

double CrossSectionDNDXIntegral::Calculate(double energy)
{
    return Calculate(energy, std::get<MAX>(GetIntegrationLimits(energy)));
}

double CrossSectionDNDXIntegral::Calculate(double energy, double v)
{
    auto lim = GetIntegrationLimits(energy);
    if (std::get<MIN>(lim) < v)
        return dndx_integral(integral, energy, std::get<MIN>(lim), v, 0);
    return 0;
}

double CrossSectionDNDXIntegral::GetUpperLimit(double energy, double rate)
{
    auto lim = GetIntegrationLimits(energy);
    auto v = dndx_upper_limit(
        integral, energy, std::get<MIN>(lim), std::get<MAX>(lim), -rate);
    return v;
}
