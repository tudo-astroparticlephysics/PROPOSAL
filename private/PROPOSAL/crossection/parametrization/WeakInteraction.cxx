
#include <cmath>

#include "PROPOSAL/Constants.h"

#include "PROPOSAL/crossection/parametrization/ParamTables.h"
#include "PROPOSAL/crossection/parametrization/WeakInteraction.h"
#include "PROPOSAL/particle/ParticleDef.h"

#include "PROPOSAL/math/Interpolant.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

WeakInteraction::WeakInteraction(
    const ParticleDef& p_def, const component_list& comp)
    : Parametrization("weak_interaction", p_def, comp, p_def.mass)
{
}

Parametrization::KinematicLimits WeakInteraction::GetKinematicLimits(
    double energy)
{
    double aux = (MP + MN) / 2; // for isoscalar targets
    aux = 2 * energy * aux + pow(aux, 2);

    auto v_min = 1e6 / (aux); // q^2_min = 1e6 MeV for experimental reasons
    auto v_max = 1;

    return KinematicLimits(v_min, v_max);
}

WeakCooperSarkarMertsch::WeakCooperSarkarMertsch(
    const ParticleDef& p_def, const component_list& comp)
    : WeakInteraction(p_def, comp)
{
    // Initialize interpolant for particles (remember crossing symmetry rules)
    if (p_def.charge < 0.) {
        interpolant_.emplace_back(new Interpolant(energies, y_nubar_p,
            sigma_nubar_p, IROMB, false, false, IROMB, false, false));
        interpolant_.emplace_back(new Interpolant(energies, y_nubar_n,
            sigma_nubar_n, IROMB, false, false, IROMB, false, false));
    } else if (p_def.charge > 0.) {
        interpolant_.emplace_back(new Interpolant(energies, y_nu_p, sigma_nu_p,
            IROMB, false, false, IROMB, false, false));
        interpolant_.emplace_back(new Interpolant(energies, y_nu_n, sigma_nu_n,
            IROMB, false, false, IROMB, false, false));
    }
}

double WeakCooperSarkarMertsch::DifferentialCrossSection(
    double energy, double v)
{
    auto log10_energy = std::log10(energy);
    const auto& nuclear_charge = current_component_.GetNucCharge();
    const auto& nuclear_number = current_component_.GetAtomicNum();

    auto proton_contr
        = nuclear_charge * interpolant_[0]->InterpolateArray(log10_energy, v);
    auto neutron_contr = (nuclear_number - nuclear_charge)
        * interpolant_[1]->InterpolateArray(log10_energy, v);
    auto mean_contr = (proton_contr + neutron_contr) / nuclear_number;

    assert(mean_contr > 0);

    // factor 1e-36: conversion from pb to cm^2
    return NA / current_component_.GetAtomicNum() * 1e-36 * mean_contr;
}
