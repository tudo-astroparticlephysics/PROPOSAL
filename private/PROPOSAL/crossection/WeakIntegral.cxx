
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/WeakIntegral.h"
#include "PROPOSAL/crossection/parametrization/WeakInteraction.h"
#include "PROPOSAL/medium/Medium.h"
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

std::pair<std::vector<Particle*>, bool> WeakIntegral::CalculateProducedParticles(double energy, double energy_loss, const Vector3D& initial_direction){
    // interaction is fatal and the initial particle is converted to a neutrino
    Particle* return_particle;
    return_particle = new Particle(*parametrization_->GetParticleDef().GetWeakPartner());
    return_particle->SetEnergy(energy - energy_loss);
    return_particle->SetDirection(initial_direction);

    return std::make_pair(std::vector<Particle*>{return_particle}, true);
}
