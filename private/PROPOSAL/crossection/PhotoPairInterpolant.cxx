#include "PROPOSAL/crossection/PhotoPairInterpolant.h"
#include "PROPOSAL/crossection/parametrization/PhotoPairProduction.h"


using namespace PROPOSAL;

PhotoPairInterpolant::PhotoPairInterpolant(unique_ptr<PhotoPairProduction>&& param,
    shared_ptr<const EnergyCutSettings> cut, const InterpolationDef& def)
    : CrossSectionInterpolant(forward<unique_ptr<PhotoPairProduction>>(param), cut, def)
{
}

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
