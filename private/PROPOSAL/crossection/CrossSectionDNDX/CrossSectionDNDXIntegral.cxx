#include "PROPOSAL/crossection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"

using namespace PROPOSAL;

double CrossSectionDNDXIntegral::Calculate(
    double energy, double v, v_trafo_t trafo)
{
    auto lim = integral_limits(energy);
    if (trafo)
        get<MAX>(lim) = trafo(get<MIN>(lim), get<MAX>(lim), v);
    return integral.Integrate(get<MIN>(lim), v,
        [&](double v) { return function_to_dndx(energy, v); }, 4);
}

double CrossSectionDNDXIntegral::GetUpperLim(
    double energy, double rnd, v_trafo_t trafo)
{
    auto lim = integral_limits(energy);
    integral.IntegrateWithRandomRatio(get<MIN>(lim), get<MAX>(lim),
        [&](double v) { return function_to_dndx(energy, v); }, 4, rnd);
    auto v = integral.GetUpperLimit();
    if (trafo)
        v = trafo(get<MIN>(lim), get<MAX>(lim), v);
    return v;
}
