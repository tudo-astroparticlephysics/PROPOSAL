
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/WeakIntegral.h"
#include "PROPOSAL/crossection/parametrization/WeakInteraction.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

WeakIntegral::WeakIntegral(const WeakInteraction& param)
        : CrossSectionIntegral(InteractionType::WeakInt, param)
{
}

WeakIntegral::WeakIntegral(const WeakIntegral& weak)
        : CrossSectionIntegral(weak)
{
}

WeakIntegral::~WeakIntegral() {}

std::pair<std::vector<DynamicData>, bool> WeakIntegral::CalculateProducedParticles(double energy, double energy_loss, const Vector3D& initial_direction){
    // interaction is fatal and the initial particle is converted to a neutrino
    DynamicData return_particle(parametrization_->GetParticleDef().weak_partner);
    // int p_id(static_cast<int>(parametrization_->GetParticleDef().weak_partner));
    // auto p_def = Id_Particle_Map.find(p_id);
    // return_particle = new Particle(p_def->second);
    return_particle.SetEnergy(energy - energy_loss);
    return_particle.SetDirection(initial_direction);

    return std::make_pair(std::vector<DynamicData>{return_particle}, true);
}
