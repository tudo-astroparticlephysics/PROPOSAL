
#include <cmath>

#include "PROPOSAL/Constants.h"

#include "PROPOSAL/crosssection/parametrization/ParamTables.h"
#include "PROPOSAL/crosssection/parametrization/WeakInteraction.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/particle/ParticleDef.h"

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

double crosssection::WeakInteraction::GetLowerEnergyLim(
    const ParticleDef& p_def) const noexcept
{
    return p_def.mass;
}

crosssection::KinematicLimits crosssection::WeakInteraction::GetKinematicLimits(
    ParticleDef const&, Component const&, double energy) const
{
    double aux = (MP + MN) / 2; // for isoscalar targets
    aux = 2 * energy * aux + pow(aux, 2);
    auto kin_lim = KinematicLimits();
    kin_lim.v_min = 1e6 / (aux); // q^2_min = 1e6 MeV for experimental reasons
    kin_lim.v_max = 1.f;
    return kin_lim;
}

std::unique_ptr<crosssection::Parametrization<Component>>
crosssection::WeakCooperSarkarMertsch::clone() const
{
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
}

crosssection::WeakCooperSarkarMertsch::WeakCooperSarkarMertsch()
{
    hash_combine(hash, std::string("cooper_sarkar_mertsch"));
    auto interpolant_particle_p = std::make_shared<Interpolant>(energies,
        y_nubar_p, sigma_nubar_p, IROMB, false, false, IROMB, false, false);
    auto interpolant_particle_n = std::make_shared<Interpolant>(energies,
        y_nubar_n, sigma_nubar_n, IROMB, false, false, IROMB, false, false);
    interpolants_particle
        = std::make_pair(interpolant_particle_p, interpolant_particle_n);

    auto interpolant_antiparticle_p = std::make_shared<Interpolant>(
        energies, y_nu_p, sigma_nu_p, IROMB, false, false, IROMB, false, false);
    auto interpolant_antiparticle_n = std::make_shared<Interpolant>(
        energies, y_nu_n, sigma_nu_n, IROMB, false, false, IROMB, false, false);
    interpolants_antiparticle = std::make_pair(
        interpolant_antiparticle_p, interpolant_antiparticle_n);
}

double crosssection::WeakCooperSarkarMertsch::DifferentialCrossSection(
    const ParticleDef& p_def, const Component& comp, double energy,
    double v) const
{
    auto log10_energy = std::log10(energy);
    const auto& nuclear_charge = comp.GetNucCharge();
    const auto& nuclear_number = comp.GetAtomicNum();

    double proton_contr = nuclear_charge;
    double neutron_contr = (nuclear_number - nuclear_charge);
    if (p_def.charge < 0.) {
        proton_contr
            *= interpolants_particle.first->InterpolateArray(log10_energy, v);
        neutron_contr
            *= interpolants_particle.second->InterpolateArray(log10_energy, v);
    } else {
        proton_contr *= interpolants_antiparticle.first->InterpolateArray(
            log10_energy, v);
        neutron_contr *= interpolants_antiparticle.second->InterpolateArray(
            log10_energy, v);
    }

    auto mean_contr = (proton_contr + neutron_contr) / nuclear_number;
    // assert(mean_contr > 0);

    // factor 1e-36: conversion from pb to cm^2
    return NA / comp.GetAtomicNum() * 1e-36 * std::max(mean_contr, 0.);
}
