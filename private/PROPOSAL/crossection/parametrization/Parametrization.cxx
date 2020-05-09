
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include <sstream>

using namespace PROPOSAL;

using std::string;
// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

Parametrization::Parametrization(const InteractionType interaction_type,
    const string& param_name, const ParticleDef& p_def,
    const component_list& components, double lower_energy_lim)
    : interaction_type_(interaction_type)
    , param_name_(param_name)
    , particle_mass_(p_def.mass)
    , particle_charge_(p_def.charge)
    , components_(components)
    , lower_energy_lim_(lower_energy_lim)
    , current_component_(components.front())
{
}

Parametrization::KinematicLimits::KinematicLimits(double v_min, double v_max)
    : vMin(v_min)
    , vMax(v_max)
{
}

double Parametrization::FunctionToDNdxIntegral(double energy, double variable)
{
    return DifferentialCrossSection(energy, variable);
}

double Parametrization::FunctionToDEdxIntegral(double energy, double variable)
{
    return variable * DifferentialCrossSection(energy, variable);
}

double Parametrization::FunctionToDE2dxIntegral(double energy, double variable)
{
    return variable * variable * DifferentialCrossSection(energy, variable);
}

size_t Parametrization::GetHash() const
{
    size_t hash_digest = 0;
    hash_combine(hash_digest, param_name_, particle_mass_, particle_charge_,
        lower_energy_lim_);
    for (const auto& comp : components_)
        hash_combine(hash_digest, comp.GetHash());

    return hash_digest;
}

void Parametrization::SetCurrentComponent(const Components::Component& comp)
{
    current_component_ = comp;
}

/* double Parametrization_builder::DifferentialCrossSection( */
/*     double energy, double v) */
/* { */
/*     return diff_cross(*this, energy, v); */
/* } */

/* Parametrization::KinematicLimits Parametrization_builder::GetKinematicLimits(
 */
/*     double energy) */
/* { */
/*     if (kinematic_limits) */
/*         return kinematic_limits(*this, energy); */

/*     return KinematicLimits(0, 0); */
/* } */

/* size_t GetHash(const std::vector<Parametrization*>& params) */
/* { */
/*     size_t hash_digest = 0; */
/*     for (auto param : params) */
/*         hash_combine(hash_digest, param->GetHash()); */
/*     return hash_digest; */
/* } */

/* namespace PROPOSAL { */
/* std::ostream& operator<<(std::ostream& os, Parametrization const& param) */
/* { */
/*     std::stringstream ss; */
/*     ss << " Parametrization (" << &param << ") "; */
/*     os << Helper::Centered(80, ss.str()) << '\n'; */

/*     os << "name: " << param.GetName() << '\n'; */

/*     param.print(os); */

/*     os << "current component index: " << param.component_index_ << '\n'; */
/*     os << "particle_mass: " << param.particle_mass_ << '\n'; */
/*     os << "particle_charge: " << param.particle_charge_ << '\n'; */
/*     os << "lower_energy_lim_: " << param.lower_energy_lim_ << '\n'; */
/*     os << Helper::Centered(80, ""); */
/*     return os; */
/* } */
/* } // namespace PROPOSAL */
