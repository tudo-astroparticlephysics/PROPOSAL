#include "PROPOSAL/crossection/PhotoPairIntegral.h"
#include "PROPOSAL/crossection/parametrization/PhotoPairProduction.h"

using namespace PROPOSAL;

PhotoPairIntegral::PhotoPairIntegral(unique_ptr<PhotoPairProduction>&& param,
    shared_ptr<const EnergyCutSettings> cut)
    : CrossSectionIntegral(forward<unique_ptr<PhotoPairProduction>>(param), cut) {}

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
