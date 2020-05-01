#include <functional>

#include "PROPOSAL/crossection/EpairInterpolant.h"
#include "PROPOSAL/crossection/parametrization/EpairProduction.h"

using namespace PROPOSAL;

EpairInterpolant::EpairInterpolant(unique_ptr<EpairProduction>&& param,
    shared_ptr<const EnergyCutSettings> cut, const InterpolationDef& interpol_def)
    : CrossSectionInterpolant(
          forward<unique_ptr<EpairProduction>>(param), cut, interpol_def)
{
}
