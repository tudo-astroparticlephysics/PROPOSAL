
#include <functional>

#include "PROPOSAL/crossection/WeakIntegral.h"
#include "PROPOSAL/crossection/WeakInterpolant.h"
#include "PROPOSAL/crossection/parametrization/WeakInteraction.h"

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

WeakInterpolant::WeakInterpolant(const WeakInteraction& param, InterpolationDef def)
        : CrossSectionInterpolant(DynamicData::WeakInt, param) {
    // Use parent CrossSecition dNdx interpolation
    InitdNdxInterpolation(def);

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

WeakInterpolant::WeakInterpolant(const WeakInterpolant& param)
        : CrossSectionInterpolant(param), converted_particle_(param.converted_particle_)
{
}

WeakInterpolant::~WeakInterpolant() {}

bool WeakInterpolant::compare(const CrossSection& cross_section) const
{
    const WeakInterpolant* cross_section_interpolant =
            static_cast<const WeakInterpolant*>(&cross_section);

    if (dndx_interpolant_1d_.size() != cross_section_interpolant->dndx_interpolant_1d_.size())
        return false;
    else if (dndx_interpolant_2d_.size() != cross_section_interpolant->dndx_interpolant_2d_.size())
        return false;

    for (unsigned int i = 0; i < dndx_interpolant_1d_.size(); ++i)
    {
        if (*dndx_interpolant_1d_[i] != *cross_section_interpolant->dndx_interpolant_1d_[i])
            return false;
    }
    for (unsigned int i = 0; i < dndx_interpolant_2d_.size(); ++i)
    {
        if (*dndx_interpolant_2d_[i] != *cross_section_interpolant->dndx_interpolant_2d_[i])
            return false;
    }

    return true;
}

std::pair<std::vector<Particle*>, bool> WeakInterpolant::CalculateProducedParticles(double energy, double energy_loss, const Vector3D initial_direction){
    // interaction is fatal and the initial particle is converted to a neutrino
    Particle* return_particle;
    return_particle = new Particle(*converted_particle_);
    return_particle->SetEnergy(energy - energy_loss);
    return_particle->SetDirection(initial_direction);

    return std::make_pair(std::vector<Particle*>{return_particle}, true);
}
