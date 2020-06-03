
#include <cmath>

#include "PROPOSAL/Constants.h"

#include "PROPOSAL/crossection/parametrization/ParamTables.h"
#include "PROPOSAL/crossection/parametrization/WeakInteraction.h"
#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/math/Interpolant.h"

using namespace PROPOSAL;
using std::make_tuple;

WeakInteraction::WeakInteraction()
    : Parametrization(InteractionType::WeakInt, "weak_interaction")
{
}

double WeakInteraction::GetLowerEnergyLim(const ParticleDef& p_def) const noexcept {
    return p_def.mass;
}
tuple<double, double> WeakInteraction::GetKinematicLimits(
    const ParticleDef& p_def, const Component& comp, double energy) const noexcept
{
    double aux = (MP + MN) / 2; // for isoscalar targets
    aux = 2 * energy * aux + pow(aux, 2);
    auto v_min = 1e6 / (aux); // q^2_min = 1e6 MeV for experimental reasons
    auto v_max = 1.f;
    return make_tuple(v_min, v_max);
}

WeakCooperSarkarMertsch::WeakCooperSarkarMertsch()
    : WeakInteraction()
{
    // Initialize interpolant for particles (remember crossing symmetry rules)
    interpolant_[false].emplace_back(new Interpolant(energies, y_nubar_p,
        sigma_nubar_p, IROMB, false, false, IROMB, false, false));
    interpolant_[false].emplace_back(new Interpolant(energies, y_nubar_n,
        sigma_nubar_n, IROMB, false, false, IROMB, false, false));
    interpolant_[true].emplace_back(new Interpolant(energies, y_nu_p,
        sigma_nu_p, IROMB, false, false, IROMB, false, false));
    interpolant_[true].emplace_back(new Interpolant(energies, y_nu_n,
        sigma_nu_n, IROMB, false, false, IROMB, false, false));
}

double WeakCooperSarkarMertsch::DifferentialCrossSection(
    const ParticleDef& p_def, const Component& comp, double energy, double v)
{
    auto log10_energy = std::log10(energy);
    const auto& nuclear_charge = comp.GetNucCharge();
    const auto& nuclear_number = comp.GetAtomicNum();

    auto proton_contr = nuclear_charge
        * interpolant_[p_def.lifetime > 0][0]->InterpolateArray(
              log10_energy, v);
    auto neutron_contr = (nuclear_number - nuclear_charge)
        * interpolant_[p_def.lifetime > 0][1]->InterpolateArray(
              log10_energy, v);
    auto mean_contr = (proton_contr + neutron_contr) / nuclear_number;

    assert(mean_contr > 0);

    // factor 1e-36: conversion from pb to cm^2
    return NA / comp.GetAtomicNum() * 1e-36 * mean_contr;
}
