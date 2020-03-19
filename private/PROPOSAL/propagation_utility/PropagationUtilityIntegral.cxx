
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

double UtilityIntegral::Calculate(double ei, double ef, double rnd)
{
    return integral.IntegrateWithRandomRatio(ei, ef,FunctionToIntegral, 4, -rnd);
}

double UtilityIntegral::GetUpperLimit(double ei, double rnd)
{
    (void)ei;
    (void)rnd;

    return integral.GetUpperLimit();
}
