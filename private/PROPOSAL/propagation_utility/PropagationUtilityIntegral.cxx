
#include <boost/bind.hpp>

#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/Constants.h"

#include "PROPOSAL/crossection/IonizIntegral.h"
#include "PROPOSAL/crossection/EpairIntegral.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"
#include "PROPOSAL/crossection/parametrization/EpairProduction.h"

using namespace PROPOSAL;

PropagationUtilityIntegral::PropagationUtilityIntegral(const ParticleDef& particle_def)
    : PropagationUtility(particle_def)
    , integral_(IROMB, IMAXS, IPREC2)
    , prop_interaction_(IROMB, IMAXS, IPREC2)
    , prop_decay_(IROMB, IMAXS, IPREC2)
    , time_particle_(IROMB, IMAXS, IPREC2)
    , dE2de_integral_(IROMB, IMAXS, IPREC2)
{
}

PropagationUtilityIntegral::PropagationUtilityIntegral(const ParticleDef& particle_def,
                                                       const Medium& medium,
                                                       const EnergyCutSettings& cut_settings,
                                                       const Definition& def)
    : PropagationUtility(particle_def, medium, cut_settings, def)
    , integral_(IROMB, IMAXS, IPREC2)
    , prop_interaction_(IROMB, IMAXS, IPREC2)
    , prop_decay_(IROMB, IMAXS, IPREC2)
    , time_particle_(IROMB, IMAXS, IPREC2)
    , dE2de_integral_(IROMB, IMAXS, IPREC2)
{
    Parametrization::Definition param_def;

    param_def.multiplier         = utility_def_.brems_multiplier;
    param_def.lpm_effect_enabled = utility_def_.lpm_effect_enabled;

    crosssections_.push_back(BremsstrahlungFactory::Get().CreateBremsstrahlung(
        utility_def_.brems_parametrization, particle_def_, *medium_, cut_settings_, param_def, false));

    param_def.multiplier = def.photo_multiplier;

    crosssections_.push_back(PhotonuclearFactory::Get().CreatePhotonuclear(
        utility_def_.photo_parametrization,
        particle_def_,
        *medium_,
        cut_settings_,
        utility_def_.photo_shadow,
        utility_def_.hardbb_enabled,
        param_def,
        false));

    param_def.multiplier = utility_def_.ioniz_multiplier;

    crosssections_.push_back(
        new IonizIntegral(Ionization(particle_def_, *medium_, cut_settings_, param_def)));

    param_def.multiplier = utility_def_.epair_multiplier;

    crosssections_.push_back(
        new EpairIntegral(EpairProductionRhoIntegral(particle_def_, *medium_, cut_settings_, param_def)));
}

PropagationUtilityIntegral::PropagationUtilityIntegral(const PropagationUtilityIntegral& collection)
    : PropagationUtility(collection)
    ,integral_(collection.integral_)
    ,prop_interaction_(collection.prop_interaction_)
    ,prop_decay_(collection.prop_decay_)
    ,time_particle_(collection.time_particle_)
    ,dE2de_integral_(collection.dE2de_integral_)
{
}

PropagationUtilityIntegral::~PropagationUtilityIntegral()
{
}

// ------------------------------------------------------------------------- //
double PropagationUtilityIntegral::CalculateDisplacement( double ei, double ef, double dist)
{
    return integral_.IntegrateWithRandomRatio(
        ei, ef, boost::bind(&PropagationUtilityIntegral::FunctionToIntegral, this,  _1), 4, -dist);
}

// ------------------------------------------------------------------------- //
double PropagationUtilityIntegral::CalculateFinalEnergy( double ei, double dist)
{
    (void)ei;
    (void)dist;

    return integral_.GetUpperLimit();
}

// ------------------------------------------------------------------------- //
double PropagationUtilityIntegral::CalculateFinalEnergy(
                                                double ei,
                                                double rnd,
                                                bool particle_interaction)
{
    (void)ei;
    (void)rnd;

    if (particle_interaction)
    {
        return prop_interaction_.GetUpperLimit();
    } else
    {
        return prop_decay_.GetUpperLimit();
    }
}

// ------------------------------------------------------------------------- //
double PropagationUtilityIntegral::CalculateTrackingIntegal(
                                                    double initial_energy,
                                                    double rnd,
                                                    bool particle_interaction)
{
    if (particle_interaction)
    {
        return prop_interaction_.IntegrateWithRandomRatio(
            initial_energy,
            particle_def_.low,
            boost::bind(&PropagationUtilityIntegral::FunctionToPropIntegralInteraction, this,  _1),
            4,
            -rnd);
    } else
    {
        return prop_decay_.IntegrateWithRandomRatio(
            initial_energy,
            particle_def_.low,
            boost::bind(&PropagationUtilityIntegral::FunctionToPropIntegralDecay, this,  _1),
            4,
            -rnd);
    }
}

// ------------------------------------------------------------------------- //
double PropagationUtilityIntegral::CalculateParticleTime( double ei, double ef)
{
    //TODO(mario): Prefactor Fri 2017/09/15
    return time_particle_.Integrate(
        ei, ef, boost::bind(&PropagationUtilityIntegral::FunctionToIntegral, this,  _1), 4);
}

// ------------------------------------------------------------------------- //
double PropagationUtilityIntegral::CalculateDE2de(double ei, double ef)
{
    return dE2de_integral_.Integrate(ei, ef, boost::bind(&PropagationUtilityIntegral::FunctionToDE2deIntegral, this, _1), 4);
}
