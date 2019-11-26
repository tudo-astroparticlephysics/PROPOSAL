#include "PROPOSAL/Secondaries.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/math/Vector3D.h"
#include <memory>
#include <vector>
#include <iostream>

using namespace PROPOSAL;

Secondaries::Secondaries() {}

Secondaries::Secondaries(size_t number_secondaries) {
    secondaries_.reserve(number_secondaries);
}

void Secondaries::push_back(const DynamicData& continuous_loss)
{
    secondaries_.push_back(continuous_loss);
}


void Secondaries::emplace_back(const InteractionType& type,const  Vector3D& position,
        const Vector3D& direction, const double& energy, const double& parent_particle_energy,
        const double& time, const double& distance)
{
    secondaries_.emplace_back(type, position, direction, energy, parent_particle_energy, time, distance);
}

void Secondaries::push_back(const Particle& particle, const InteractionType& interaction_type, const double& energy_loss)
{
    DynamicData data(interaction_type);

    data.SetEnergy(energy_loss);
    data.SetPosition(particle.GetPosition());
    data.SetDirection(particle.GetDirection());
    data.SetTime(particle.GetTime());
    data.SetParentParticleEnergy(particle.GetEnergy());
    data.SetPropagatedDistance(particle.GetPropagatedDistance());

    secondaries_.push_back(data);
}

void Secondaries::append(Secondaries secondaries)
{
    secondaries_.shrink_to_fit();
    secondaries_.insert(secondaries_.end(), secondaries.secondaries_.begin(), secondaries.secondaries_.end());
}
