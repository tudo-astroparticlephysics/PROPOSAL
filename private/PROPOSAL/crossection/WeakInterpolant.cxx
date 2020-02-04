
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
        : CrossSectionInterpolant(InteractionType::WeakInt, param) {
    // Use parent CrossSecition dNdx interpolation
    InitdNdxInterpolation(def);
}

WeakInterpolant::WeakInterpolant(const WeakInterpolant& param)
        : CrossSectionInterpolant(param)
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

std::pair<std::vector<DynamicData>, bool> WeakInterpolant::CalculateProducedParticles(double energy, double energy_loss, const Vector3D& initial_direction){
    // interaction is fatal and the initial particle is converted to a neutrino
    DynamicData return_particle(parametrization_->GetParticleDef().weak_partner);
    // int p_id(static_cast<int>(parametrization_->GetParticleDef().weak_partner));
    // auto p_def = Id_Particle_Map.find(p_id);
    // return_particle = new Particle(p_def->second);
    return_particle.SetEnergy(energy - energy_loss);
    return_particle.SetDirection(initial_direction);

    return std::make_pair(std::vector<DynamicData>{return_particle}, true);
}
