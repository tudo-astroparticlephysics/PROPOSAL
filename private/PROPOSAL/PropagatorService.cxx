
#include "PROPOSAL/PropagatorService.h"
#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/Secondaries.h"

#include <memory>

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
PropagatorService::PropagatorService()
    : propagator_map_()
{
}

// ------------------------------------------------------------------------- //
PropagatorService::~PropagatorService()
{
    for (PropagatorMap::const_iterator it = propagator_map_.begin(); it != propagator_map_.end(); ++it)
    {
        delete it->second;
    }

    propagator_map_.clear();
}

// ------------------------------------------------------------------------- //
void PropagatorService::RegisterPropagator(const Propagator& propagator)
{
    ParticleDef particle_def = const_cast<Propagator&>(propagator).GetParticleDef();
    if (propagator_map_.find(particle_def) != propagator_map_.end())
    {
        log_warn("Propagator for particle %s is already registered!", particle_def.name.c_str());
    } else
    {
        propagator_map_[particle_def] = new Propagator(propagator);
    }
}

// ------------------------------------------------------------------------- //
bool PropagatorService::IsRegistered(const ParticleDef& particle_def)
{
    return propagator_map_.find(particle_def) != propagator_map_.end();
}

// ------------------------------------------------------------------------- //
Secondaries PropagatorService::Propagate(
    const ParticleDef& particle_def,
    DynamicData& particle_condition,
    double max_distance,
    double min_energy)
{
    PropagatorMap::iterator it = propagator_map_.find(particle_def);

    if (it != propagator_map_.end())
    {
        Propagator* propagator = it->second;

        // Particle& prop_particle = propagator->GetParticle();
        // prop_particle.SetEnergy(particle.GetEnergy());
        // prop_particle.SetParentParticleEnergy(particle.GetParentParticleEnergy());
        // prop_particle.SetPosition(particle.GetPosition());
        // prop_particle.SetDirection(particle.GetDirection());
        // prop_particle.SetPropagatedDistance(particle.GetPropagatedDistance());
        // prop_particle.SetTime(particle.GetTime());
        // prop_particle.SetElost(particle.GetElost());

        return propagator->Propagate(particle_condition, max_distance, min_energy);

        // particle.SetParentParticleEnergy(prop_particle.GetParentParticleEnergy());
        // particle_condition.SetPropagatedDistance(prop_particle.GetPropagatedDistance());
        // particle.SetElost(prop_particle.GetElost());
        // particle.SetEntryPoint(prop_particle.GetEntryPoint());
        // particle.SetEntryEnergy(prop_particle.GetEntryEnergy());
        // particle.SetEntryTime(prop_particle.GetEntryTime());
        // particle.SetExitPoint(prop_particle.GetExitPoint());
        // particle.SetExitEnergy(prop_particle.GetExitEnergy());
        // particle.SetExitTime(prop_particle.GetExitTime());
        // particle.SetClosestApproachPoint(prop_particle.GetClosestApproachPoint());
        // particle.SetClosestApproachEnergy(prop_particle.GetClosestApproachEnergy());
        // particle.SetClosestApproachTime(prop_particle.GetClosestApproachTime());
        // particle.InjectState(prop_particle);

    } else
    {
        log_warn("Propagator for particle %s not found! Empty secondary vector will be returned!",
                 particle_def.name.c_str());
        return Secondaries();
    }
}

Propagator* PropagatorService::GetPropagatorToParticleDef(const ParticleDef& particle_def)
{
    PropagatorMap::iterator it = propagator_map_.find(particle_def);
    if (it != propagator_map_.end())
    {
        return it->second;
    }
    else
    {
        log_warn("Propagator for particle %s not found! Default nullptr will be returned!",
                 particle_def.name.c_str());
        Propagator* default_propagator = nullptr;
        return default_propagator;
    }
}
