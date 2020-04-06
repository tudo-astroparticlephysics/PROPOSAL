
#include <cmath>
#include <functional>

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

UtilityIntegral::UtilityIntegral(std::function<double(double)> func, double lower_lim)
    : integral(IROMB, IMAXS, IPREC2)
    , FunctionToIntegral(func)
    , lower_lim(lower_lim)
{
}

void UtilityIntegral::BuildTables(const std::string str, size_t hash_digest,
    Interpolant1DBuilder::Definition interpol_def)
{
    (void)str;
    (void)hash_digest;
    (void)interpol_def;
};

double UtilityIntegral::Calculate(
    double energy_initial, double energy_final, double rnd)
{
    assert(energy_initial > energy_final);
    assert(energy_final > lower_lim);

    last_energy_initial = energy_initial;
    last_partial_sum = rnd;

    return integral.IntegrateWithRandomRatio(
        energy_initial, energy_final, FunctionToIntegral, 4, -rnd);
}

double UtilityIntegral::GetUpperLimit(double energy_initial, double rnd)
{
    if (energy_initial != last_energy_initial || rnd != last_partial_sum)
        Calculate(energy_initial, lower_lim, rnd);

    return integral.GetUpperLimit();
}
