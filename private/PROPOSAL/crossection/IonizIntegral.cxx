
#include <functional>

#include <cmath>

#include "PROPOSAL/crossection/IonizIntegral.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"

#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;


double IonizIntegral::dedx_integral(double energy)
{
    if (parametrization_->GetName() == "IonizBergerSeltzerBhabha"
        || parametrization_->GetName() == "IonizBergerSeltzerMoller") {
        return parametrization_->FunctionToDEdxIntegral(energy, 0);
    }

    return CrossSectionIntegral::dedx_integral(energy);
}
