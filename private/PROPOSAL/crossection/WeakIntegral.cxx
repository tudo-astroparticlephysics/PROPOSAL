
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/WeakIntegral.h"
#include "PROPOSAL/crossection/parametrization/WeakInteraction.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

WeakIntegral::WeakIntegral(const WeakInteraction& param)
        : CrossSectionIntegral(DynamicData::WeakInt, param)
{
    const ParticleDef& param_particledef = parametrization_->GetParticleDef();

    if(param_particledef==EMinusDef::Get()){
        converted_particle_ = &NuEDef::Get();
    }
    else if (param_particledef==MuMinusDef::Get()){
        converted_particle_ = &NuMuDef::Get();
    }
    else if (param_particledef==TauMinusDef::Get()){
        converted_particle_ = &NuTauDef::Get();
    }
    else if (param_particledef==EPlusDef::Get()){
        converted_particle_ = &NuEBarDef::Get();
    }
    else if (param_particledef==MuPlusDef::Get()){
        converted_particle_ = &NuMuBarDef::Get();
    }
    else if (param_particledef==TauPlusDef::Get()){
        converted_particle_ = &NuTauBarDef::Get();
    }
    else{
        log_fatal("Weak interaction: Particle to propagate is not a SM charged lepton");
    }
}

WeakIntegral::WeakIntegral(const WeakIntegral& weak)
        : CrossSectionIntegral(weak), converted_particle_(weak.converted_particle_)
{
}

WeakIntegral::~WeakIntegral() {}

std::pair<std::vector<Particle*>, bool> WeakIntegral::CalculateProducedParticles(double energy, double energy_loss, const Vector3D initial_direction){
    // interaction is fatal and the initial particle is converted to a neutrino
    Particle* return_particle;
    return_particle = new Particle(*converted_particle_);
    return_particle->SetEnergy(energy - energy_loss);
    return_particle->SetDirection(initial_direction);

    return std::make_pair(std::vector<Particle*>{return_particle}, true);
}