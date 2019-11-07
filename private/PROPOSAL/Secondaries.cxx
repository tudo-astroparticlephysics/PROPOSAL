#include "PROPOSAL/Secondaries.h"
#include "PROPOSAL/particle/Particle.h"
#include <memory>

using namespace PROPOSAL;

void Secondaries::push_back(const Particle& particle)
{
    Particle data(particle);
    secondarys_.push_back(std::make_shared<DynamicData>(data));
}

void Secondaries::push_back(DynamicData continuous_loss)
{
    secondarys_.push_back(std::make_shared<DynamicData>(continuous_loss));
}

void Secondaries::push_back(std::shared_ptr<DynamicData> continuous_loss){
    secondarys_.push_back(continuous_loss);
}

void Secondaries::push_back(const Particle& particle, const DynamicData::Type& interaction_type, double energy_loss)
{
    DynamicData data(interaction_type);

    data.SetEnergy(energy_loss);
    data.SetPosition(particle.GetPosition());
    data.SetDirection(particle.GetDirection());
    data.SetTime(particle.GetTime());
    data.SetParentParticleEnergy(particle.GetEnergy());
    data.SetPropagatedDistance(particle.GetPropagatedDistance());

    secondarys_.push_back(std::make_shared<DynamicData>(data));
}
