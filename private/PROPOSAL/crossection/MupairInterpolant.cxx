#include "PROPOSAL/crossection/MupairInterpolant.h"
#include "PROPOSAL/crossection/parametrization/MupairProduction.h"

using namespace PROPOSAL;

MupairInterpolant::MupairInterpolant(unique_ptr<MupairProduction>&& param,
    shared_ptr<const EnergyCutSettings> cut,
    const InterpolationDef& interpol_def)
    : CrossSectionInterpolant(forward<unique_ptr<MupairProduction>>(param), cut, interpol_def)
{
}
