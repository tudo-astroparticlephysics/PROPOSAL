#include <functional>

#include "PROPOSAL/crossection/BremsInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"

using namespace PROPOSAL;

BremsInterpolant::BremsInterpolant(unique_ptr<Bremsstrahlung>&& param,
    shared_ptr<const EnergyCutSettings> cuts,
    const InterpolationDef& interpol_def)
    : CrossSectionInterpolant(forward<unique_ptr<Bremsstrahlung>>(param), cuts, interpol_def)
{
}
