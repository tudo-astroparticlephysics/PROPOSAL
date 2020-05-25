
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


/* double IonizIntegral::CalculatedEdx( */
/*     const ParticleDef& p_def, const Medium& medium, double energy) */
/* { */
/*     auto is_bhabha = parametrization_->name == "IonizBergerSeltzerBhabha"; */
/*     auto is_moller = parametrization_->name == "IonizBergerSeltzerMoller"; */
/*     if (is_bhabha || is_moller) */
/*         return parametrization_->FunctionToDEdxIntegral( */
/*             p_def, medium, energy, 0.); */

/*     auto physical_lim = reinterpret_cast<Ionization*>(parametrization_.get()) */
/*                             ->GetKinematicLimits(p_def, medium, energy); */
/*     auto v_cut = GetEnergyCut(energy, physical_lim); */
/*     auto dEdx = [this, &p_def, &medium, energy](double v) { */
/*         return parametrization_->FunctionToDEdxIntegral( */
/*             p_def, medium, energy, v); */
/*     }; */
/*     return integral_.Integrate(physical_lim.vMin, v_cut, dEdx, 4); */
/* } */
