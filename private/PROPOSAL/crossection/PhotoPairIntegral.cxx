
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/PhotoPairIntegral.h"
#include "PROPOSAL/crossection/parametrization/PhotoPairProduction.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Logging.h"


using namespace PROPOSAL;

double PhotoPairIntegral::CalculatedEdx(double energy)
{
        (void)energy;
        return 0;
}

double PhotoPairIntegral::CalculatedE2dx(double energy)
{
        (void)energy;
        return 0;
}
