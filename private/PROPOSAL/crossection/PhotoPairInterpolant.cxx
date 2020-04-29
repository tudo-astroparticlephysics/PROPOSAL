
#include <cmath>
#include <functional>

#include "PROPOSAL/crossection/PhotoPairIntegral.h"
#include "PROPOSAL/crossection/PhotoPairInterpolant.h"
#include "PROPOSAL/crossection/parametrization/PhotoPairProduction.h"

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;
double PhotoPairInterpolant::CalculatedEdx(double energy)
{
    (void)energy;
    return 0;
}

double PhotoPairInterpolant::CalculatedE2dx(double energy)
{
    (void)energy;
    return 0;
}
