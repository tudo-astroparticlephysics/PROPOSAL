
#include "PROPOSAL/PropagatorService.h"
#include "PROPOSAL/Propagator.h"

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
    ParticleDef particle_def = const_cast<Propagator&>(propagator).GetParticle().GetParticleDef();
    if (propagator_map_.find(particle_def) != propagator_map_.end())
    {
        log_warn("Propagator for particle %s is already registered!", particle_def.name.c_str());
    }
    else
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
std::vector<DynamicData*> PropagatorService::Propagate(Particle& particle)
{
    ParticleDef particle_def = particle.GetParticleDef();

    PropagatorMap::iterator it = propagator_map_.find(particle_def);

    if (it != propagator_map_.end())
    {
        Propagator* propagator = it->second;

        Particle& prop_particle = propagator->GetParticle();
        prop_particle.SetEnergy(particle.GetEnergy());
        prop_particle.SetParentParticleEnergy(particle.GetParentParticleEnergy());
        prop_particle.SetPosition(particle.GetPosition());
        prop_particle.SetDirection(particle.GetDirection());
        prop_particle.SetPropagatedDistance(particle.GetPropagatedDistance());
        prop_particle.SetTime(particle.GetTime());
        prop_particle.SetElost(particle.GetElost());

        std::vector<DynamicData*> secondaries = propagator->Propagate();

        particle.SetEnergy(prop_particle.GetEnergy());
        particle.SetParentParticleEnergy(prop_particle.GetParentParticleEnergy());
        particle.SetPosition(prop_particle.GetPosition());
        particle.SetDirection(prop_particle.GetDirection());
        particle.SetPropagatedDistance(prop_particle.GetPropagatedDistance());
        particle.SetTime(prop_particle.GetTime());
        particle.SetElost(prop_particle.GetElost());

        return secondaries;
    }
    else
    {
        log_warn("Propagator for particle %s not found! Empty secondary vector will be returned!", particle_def.name.c_str());
        return std::vector<DynamicData*>();
    }
}
