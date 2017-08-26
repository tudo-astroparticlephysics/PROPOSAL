
#include <boost/bind.hpp>

#include "PROPOSAL/CollectionIntegral.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

CollectionIntegral::CollectionIntegral()
    : Collection()
    , integral_(IROMB, IMAXS, IPREC2)
    , prop_interaction_(IROMB, IMAXS, IPREC2)
    , prop_decay_(IROMB, IMAXS, IPREC2)
    , time_particle_(IROMB, IMAXS, IPREC2)
{
}

CollectionIntegral::CollectionIntegral(const Medium& medium,
                                       const Geometry& geometry,
                                       const EnergyCutSettings& cut_settings,
                                       const CollectionDef& def)
    : Collection(medium, geometry, cut_settings, def)
    , integral_(IROMB, IMAXS, IPREC2)
    , prop_interaction_(IROMB, IMAXS, IPREC2)
    , prop_decay_(IROMB, IMAXS, IPREC2)
    , time_particle_(IROMB, IMAXS, IPREC2)
{
}

CollectionIntegral::CollectionIntegral(const CollectionIntegral& collection)
    : Collection(collection)
    ,integral_(collection.integral_)
    ,prop_interaction_(collection.prop_interaction_)
    ,prop_decay_(collection.prop_decay_)
    ,time_particle_(collection.time_particle_)
{
}

CollectionIntegral::~CollectionIntegral()
{
}

// ------------------------------------------------------------------------- //
double CollectionIntegral::CalculateDisplacement(const PROPOSALParticle& particle, double ei, double ef, double dist)
{
    return integral_.IntegrateWithRandomRatio(
        ei, ef, boost::bind(&CollectionIntegral::FunctionToIntegral, this, boost::cref(particle), _1), 4, -dist);
}

// ------------------------------------------------------------------------- //
double CollectionIntegral::CalculateFinalEnergy(const PROPOSALParticle& particle, double ei, double dist)
{
    (void)particle;
    (void)ei;
    (void)dist;

    return integral_.GetUpperLimit();
}

// ------------------------------------------------------------------------- //
double CollectionIntegral::CalculateFinalEnergy(const PROPOSALParticle& particle,
                                                double ei,
                                                double rnd,
                                                bool particle_interaction)
{
    (void)particle;
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
double CollectionIntegral::CalculateTrackingIntegal(const PROPOSALParticle& particle,
                                                    double initial_energy,
                                                    double rnd,
                                                    bool particle_interaction)
{
    if (particle_interaction)
    {
        return prop_interaction_.IntegrateWithRandomRatio(
            initial_energy,
            particle.GetLow(),
            boost::bind(&CollectionIntegral::FunctionToPropIntegralInteraction, this, boost::cref(particle), _1),
            4,
            -rnd);
    } else
    {
        return prop_decay_.IntegrateWithRandomRatio(
            initial_energy,
            particle.GetLow(),
            boost::bind(&CollectionIntegral::FunctionToPropIntegralDecay, this, boost::cref(particle), _1),
            4,
            -rnd);
    }
}

// ------------------------------------------------------------------------- //
double CollectionIntegral::CalculateParticleTime(const PROPOSALParticle& particle, double ei, double ef)
{
    return time_particle_.Integrate(
        ei, ef, boost::bind(&CollectionIntegral::FunctionToTimeIntegral, this, boost::cref(particle), _1), 4);
}
