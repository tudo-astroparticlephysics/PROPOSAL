
#include <cmath>

#include "PROPOSAL/crossection/parametrization/WeakInteraction.h"

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/crossection/parametrization/WeakTable.h"


using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

WeakInteraction::WeakInteraction(const ParticleDef& particle_def,
                                 const Medium& medium,
                                 double multiplier)
        : Parametrization(particle_def, medium, EnergyCutSettings(), multiplier)
{
}

WeakInteraction::WeakInteraction(const WeakInteraction& param)
        : Parametrization(param)
{
}

WeakInteraction::~WeakInteraction() {}

bool WeakInteraction::compare(const Parametrization& parametrization) const
{
    return Parametrization::compare(parametrization);
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

Parametrization::IntegralLimits WeakInteraction::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    double aux = (MP + MN) / 2; //for isoscalar targets
    aux = 2 * energy * aux + pow(aux, 2);

    limits.vMin = 1e6 / (aux ); //q^2_min = 1e6 MeV for experimental reasons
    limits.vUp = limits.vMin; //treat interaction as fully stochastic
    limits.vMax = 1;

    return limits;
}

// ------------------------------------------------------------------------- //
// Specific implementations
// ------------------------------------------------------------------------- //

WeakCooperSarkarMertsch::WeakCooperSarkarMertsch(const ParticleDef& particle_def,
                                                 const Medium& medium,
                                                 double multiplier)
        : WeakInteraction(particle_def, medium, multiplier)
        , interpolant_(2, NULL)
{

    std::vector<Interpolant2DBuilder_array_as> builder2d(2); // [0] for proton, [1] for neutron
    Helper::InterpolantBuilderContainer builder_container2d(2);
    if(particle_def==EMinusDef::Get()||particle_def==MuMinusDef::Get()||particle_def==TauMinusDef::Get()){
        builder2d[0]
                .Setx1(energies)
                .Setx2(y_nubar_p)
                .Sety(sigma_nubar_p)
                .SetRomberg1(IROMB)
                .SetRational1(false)
                .SetRelative1(false)
                .SetRomberg2(IROMB)
                .SetRational2(false)
                .SetRelative2(false);

        builder_container2d[0].first  = &builder2d[0];
        builder_container2d[0].second = &interpolant_[0];

        builder2d[1]
                .Setx1(energies)
                .Setx2(y_nubar_n)
                .Sety(sigma_nubar_n)
                .SetRomberg1(IROMB)
                .SetRational1(false)
                .SetRelative1(false)
                .SetRomberg2(IROMB)
                .SetRational2(false)
                .SetRelative2(false);

        builder_container2d[1].first  = &builder2d[1];
        builder_container2d[1].second = &interpolant_[1];
    }
    else if(particle_def==EPlusDef::Get()||particle_def==MuPlusDef::Get()||particle_def==TauPlusDef::Get()){
        builder2d[0]
                .Setx1(energies)
                .Setx2(y_nu_p)
                .Sety(sigma_nu_p)
                .SetRomberg1(IROMB)
                .SetRational1(false)
                .SetRelative1(false)
                .SetRomberg2(IROMB)
                .SetRational2(false)
                .SetRelative2(false);

        builder_container2d[0].first  = &builder2d[0];
        builder_container2d[0].second = &interpolant_[0];

        builder2d[1]
                .Setx1(energies)
                .Setx2(y_nu_n)
                .Sety(sigma_nu_n)
                .SetRomberg1(IROMB)
                .SetRational1(false)
                .SetRelative1(false)
                .SetRomberg2(IROMB)
                .SetRational2(false)
                .SetRelative2(false);

        builder_container2d[1].first  = &builder2d[1];
        builder_container2d[1].second = &interpolant_[1];
    }
    else{
        log_fatal("Weak interaction table_read: Particle to propagate is not a SM charged lepton");
    }
    Helper::InitializeInterpolation("WeakInt", builder_container2d, std::vector<Parametrization*>(1, this), InterpolationDef());

}

WeakCooperSarkarMertsch::WeakCooperSarkarMertsch(const WeakCooperSarkarMertsch& param)
        : WeakInteraction(param)
        , interpolant_()
{
    interpolant_.resize(param.interpolant_.size());

    for (unsigned int i = 0; i < param.interpolant_.size(); ++i)
    {
        interpolant_[i] = new Interpolant(*param.interpolant_[i]);
    }
}

WeakCooperSarkarMertsch::~WeakCooperSarkarMertsch()
{
    for (std::vector<Interpolant*>::const_iterator iter = interpolant_.begin(); iter != interpolant_.end(); ++iter)
    {
        delete *iter;
    }

    interpolant_.clear();
}

bool WeakCooperSarkarMertsch::compare(const Parametrization& parametrization) const
{
    const WeakCooperSarkarMertsch* weak = static_cast<const WeakCooperSarkarMertsch*>(&parametrization);

    if (interpolant_.size() != weak->interpolant_.size())
        return false;

    for (unsigned int i = 0; i < interpolant_.size(); ++i)
    {
        if (*interpolant_[i] != *weak->interpolant_[i])
            return false;
    }

    return WeakInteraction::compare(parametrization);
}

double WeakCooperSarkarMertsch::DifferentialCrossSection(double energy, double v)
{
    double proton_contribution = components_[component_index_]->GetNucCharge() * interpolant_.at(0)->InterpolateArray(std::log10(energy), v);
    double neutron_contribution =  (components_[component_index_]->GetAtomicNum() - components_[component_index_]->GetNucCharge()) * interpolant_.at(1)->InterpolateArray(std::log10(energy), v);
    double mean_contribution = (proton_contribution + neutron_contribution) / (components_[component_index_]->GetAtomicNum());

    return medium_->GetMolDensity() * components_[component_index_]->GetAtomInMolecule() * 1e-36 * std::max(0.0, mean_contribution); //factor 1e-36: conversion from pb to cm^2
}


const std::string WeakCooperSarkarMertsch::name_ = "WeakCooperSarkarMertsch";

