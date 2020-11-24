
#include <cmath>
#include <functional>

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

UtilityIntegral::UtilityIntegral(
    std::function<double(double)> func, double lower_lim)
    : integral(IROMB, IMAXS, IPREC2)
    , lower_lim(lower_lim)
    , FunctionToIntegral(func)
{
}

void UtilityIntegral::BuildTables(const std::string , size_t ,
    Interpolant1DBuilder::Definition, bool)
{
}

double UtilityIntegral::Calculate(double energy_initial, double energy_final)
{
    return integral.Integrate(
        energy_initial, energy_final, FunctionToIntegral, 4);
}

double UtilityIntegral::GetUpperLimit(double energy_initial, double rnd)
{
    auto sum = integral.IntegrateWithRandomRatio(
        energy_initial, lower_lim, FunctionToIntegral, 4, -rnd);

    assert(sum > rnd); // searched Energy is below lower_lim return lower_lim as
                       // a lower limit
    (void)sum;

    return integral.GetUpperLimit();
}
