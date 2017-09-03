
#include <boost/bind.hpp>

#include "PROPOSAL/sector/SectorIntegral.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

SectorIntegral::SectorIntegral(PROPOSALParticle& particle)
    : Sector(particle)
    , integral_(IROMB, IMAXS, IPREC2)
    , prop_interaction_(IROMB, IMAXS, IPREC2)
    , prop_decay_(IROMB, IMAXS, IPREC2)
    , time_particle_(IROMB, IMAXS, IPREC2)
{
}

SectorIntegral::SectorIntegral(PROPOSALParticle& particle, const Medium& medium,
                                       const Geometry& geometry,
                                       const EnergyCutSettings& cut_settings,
                                       const Definition& def)
    : Sector(particle, medium, geometry, cut_settings, def)
    , integral_(IROMB, IMAXS, IPREC2)
    , prop_interaction_(IROMB, IMAXS, IPREC2)
    , prop_decay_(IROMB, IMAXS, IPREC2)
    , time_particle_(IROMB, IMAXS, IPREC2)
{
}

SectorIntegral::SectorIntegral(const SectorIntegral& collection)
    : Sector(collection)
    ,integral_(collection.integral_)
    ,prop_interaction_(collection.prop_interaction_)
    ,prop_decay_(collection.prop_decay_)
    ,time_particle_(collection.time_particle_)
{
}

SectorIntegral::~SectorIntegral()
{
}

// ------------------------------------------------------------------------- //
double SectorIntegral::CalculateDisplacement( double ei, double ef, double dist)
{
    return integral_.IntegrateWithRandomRatio(
        ei, ef, boost::bind(&SectorIntegral::FunctionToIntegral, this,  _1), 4, -dist);
}

// ------------------------------------------------------------------------- //
double SectorIntegral::CalculateFinalEnergy( double ei, double dist)
{
    (void)ei;
    (void)dist;

    return integral_.GetUpperLimit();
}

// ------------------------------------------------------------------------- //
double SectorIntegral::CalculateFinalEnergy(
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
double SectorIntegral::CalculateTrackingIntegal(
                                                    double initial_energy,
                                                    double rnd,
                                                    bool particle_interaction)
{
    if (particle_interaction)
    {
        return prop_interaction_.IntegrateWithRandomRatio(
            initial_energy,
            particle_.GetLow(),
            boost::bind(&SectorIntegral::FunctionToPropIntegralInteraction, this,  _1),
            4,
            -rnd);
    } else
    {
        return prop_decay_.IntegrateWithRandomRatio(
            initial_energy,
            particle_.GetLow(),
            boost::bind(&SectorIntegral::FunctionToPropIntegralDecay, this,  _1),
            4,
            -rnd);
    }
}

// ------------------------------------------------------------------------- //
double SectorIntegral::CalculateParticleTime( double ei, double ef)
{
    return time_particle_.Integrate(
        ei, ef, boost::bind(&SectorIntegral::FunctionToTimeIntegral, this,  _1), 4);
}
