
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/EpairIntegral.h"
#include "PROPOSAL/crossection/parametrization/EpairProduction.h"
#include "PROPOSAL/medium/Medium.h"

using namespace PROPOSAL;

/* EpairIntegral::EpairIntegral(unique_ptr<EpairProduction>&& param, */
/*     shared_ptr<const EnergyCutSettings> cut) */
/*     : CrossSectionIntegral(forward<unique_ptr<EpairProduction>>(param), cut) */
/* { */
/* } */

/* double EpairIntegral::FunctionToDEdxIntegralReverse( */
/*     const ParticleDef& p_def, const Component& comp, double energy, double v) */
/* { */
/*     return (1 - v) * parametrization_->DifferentialCrossSection(p_def, comp, energy, v); */
/* } */

/* double EpairIntegral::integrate_dedx( */
/*     const ParticleDef& p_def, const Component& comp, double energy) */
/* { */
/*     auto physical_lim= parametrization_->GetKinematicLimits(p_def, comp, energy); */
/*     auto v_cut = GetEnergyCut(energy, physical_lim); */

/*     double r1 = 0.8; */
/*     double rUp = v_cut * (1 - HALF_PRECISION); */
/*     bool rflag = false; */

/*     if (r1 < rUp) { */
/*         if (2 * parametrization_->FunctionToDEdxIntegral(p_def, comp, energy, r1) */
/*             < parametrization_->FunctionToDEdxIntegral(p_def, comp, energy, rUp)) { */
/*             rflag = true; */
/*         } */
/*     } */

/*     auto func = [&, energy](double v) { */
/*         return parametrization_->FunctionToDEdxIntegral(p_def, comp, energy, v); */
/*     }; */

/*     if (rflag) { */
/*         if (r1 > v_cut) { */
/*             r1 = v_cut; */
/*         } */

/*         if (r1 < physical_lim.vMin) { */
/*             r1 = physical_lim.vMin; */
/*         } */

/*         auto sum = integral_.Integrate(physical_lim.vMin, r1, func, 4); */
/*         double r2 = std::max(1 - v_cut, COMPUTER_PRECISION); */

/*         if (r2 > 1 - r1) { */
/*             r2 = 1 - r1; */
/*         } */

/*         auto func_reverse = [&, energy](double v) { */
/*             return FunctionToDEdxIntegralReverse(p_def, comp, energy, v); */
/*         }; */
/*         sum += integral_.Integrate(1 - v_cut, r2, func_reverse, 2) */
/*             + integral_.Integrate(r2, 1 - r1, func_reverse, 4); */
/*         return energy * sum; */
/*     } */
/*     return energy * integral_.Integrate(physical_lim.vMin, v_cut, func, 4); */
/* } */
