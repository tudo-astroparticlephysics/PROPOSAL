
#include <functional>

#include "PROPOSAL/crossection/IonizInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"

using namespace PROPOSAL;

IonizInterpolant::IonizInterpolant(unique_ptr<Ionization>&& param, shared_ptr<const EnergyCutSettings> cuts, const InterpolationDef& interpol_def)
    : CrossSectionInterpolant(forward<unique_ptr<Ionization>>(param), cuts, interpol_def)
{
}
