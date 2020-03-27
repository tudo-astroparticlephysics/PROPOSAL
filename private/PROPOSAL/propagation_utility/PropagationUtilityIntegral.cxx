
#include <cmath>
#include <functional>

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

UtilityIntegral::UtilityIntegral(std::function<double(double)> func)
    : integral(IROMB, IMAXS, IPREC2), FunctionToIntegral(func)
{
}

double UtilityIntegral::Calculate(double energy_initial, double energy_final, double rnd)
{
    last_energy_initial = energy_initial;
    last_partial_sum = rnd;

    return integral.IntegrateWithRandomRatio(energy_initial, energy_final, FunctionToIntegral, 4, -rnd);
}

double UtilityIntegral::GetUpperLimit(double energy_initial, double rnd)
{
    if (energy_initial != last_energy_initial || rnd != last_partial_sum)
        Calculate(energy_initial , 1e-3, rnd);

    return integral.GetUpperLimit();
}
