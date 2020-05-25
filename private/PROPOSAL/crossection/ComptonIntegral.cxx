
#include <cmath>
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/ComptonIntegral.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"

using namespace PROPOSAL;
using std::exp;
using std::log;
using std::vector;

/* ComptonIntegral::ComptonIntegral( */
/*     unique_ptr<Compton>&& param, shared_ptr<const EnergyCutSettings> cut) */
/*     : CrossSectionIntegral(forward<unique_ptr<Compton>>(param), cut) */
/* { */
/* } */

/* // Integrate with the substitution t = ln(1-v) to avoid numerical problems */
/* double ComptonIntegral::log_substitution(double v) const { return log(1 - v); } */

/* double ComptonIntegral::dndx_integral(double energy, double rnd) */
/* { */
/*     auto v_cut = GetEnergyCut(energy); */
/*     auto v_max = parametrization_->GetKinematicLimits(energy).vMax; */

/*     auto t_min = log_substitution(v_cut); */
/*     auto t_max = log_substitution(v_max); */

/*     /1* auto dndx_func = bind(&Compton::FunctionToDNdxIntegral, parametrization_, */
/*      * energy, _1); *1/ */
/*     auto func_transformed = [&](double t) { */
/*         return exp(t) */
/*             * parametrization_->FunctionToDNdxIntegral(energy, 1 - exp(t)); */
/*     }; */

/*     integral_.IntegrateWithRandomRatio(t_min, t_max, func_transformed, 4, rnd); */

/*     return integral_.GetUpperLimit(); */
/* } */

/* double ComptonIntegral::dedx_integral(double energy) */
/* { */
/*     auto v_cut = GetEnergyCut(energy); */
/*     auto v_max = parametrization_->GetKinematicLimits(energy).vMax; */

/*     auto t_min = log_substitution(v_cut); */
/*     auto t_max = log_substitution(v_max); */

/*     /1* auto dedx_func = bind( *1/ */
/*     /1*     &Parametrization::FunctionToDEdxIntegral, parametrization_, energy, */
/*      * _1); *1/ */
/*     auto func_transformed = [&](double t) { */
/*         return exp(t) */
/*             * (parametrization_->FunctionToDEdxIntegral)(energy, 1 - exp(t)); */
/*     }; */

/*     return integral_.Integrate(t_min, t_max, func_transformed, 2); */
/* } */

/* double ComptonIntegral::de2dx_integral(double energy) */
/* { */
/*     auto v_cut = GetEnergyCut(energy); */
/*     auto v_max = parametrization_->GetKinematicLimits(energy).vMax; */

/*     auto t_min = log_substitution(v_cut); */
/*     auto t_max = log_substitution(v_max); */

/*     /1* auto de2dx_func = bind(&Parametrization::FunctionToDE2dxIntegral, *1/ */
/*     /1*     parametrization_, energy, _1); *1/ */
/*     auto func_transformed = [&](double t) { */
/*         return exp(t) */
/*             * (parametrization_->FunctionToDE2dxIntegral)(energy, 1 - exp(t)); */
/*     }; */

/*     return integral_.Integrate(t_min, t_max, func_transformed, 2); */
/* } */
