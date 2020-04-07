
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include <cmath>
#include <sstream>

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

Parametrization::Parametrization(const ParticleDef& particle_def,
    std::shared_ptr<const Medium> medium, double lower_energy_lim, double multiplier)
    : particle_mass_(particle_def.mass)
    , particle_charge_(particle_def.charge)
    , particle_lifetime_(particle_def.lifetime)
    , lower_energy_lim_(lower_energy_lim)
    , medium_(medium)
    , components_(medium_->GetComponents())
    , component_index_(0)
    , multiplier_(multiplier)
{
}

Parametrization::Parametrization(const Parametrization& param)
    : particle_mass_(param.particle_mass_)
    , particle_charge_(param.particle_charge_)
    , lower_energy_lim_(param.lower_energy_lim_)
    , medium_(param.medium_)
    , components_(medium_->GetComponents())
    , component_index_(param.component_index_) // //TODO(mario): Check better
                                               // way Mon 2017/09/04
    , multiplier_(param.multiplier_)
{
}

Parametrization::~Parametrization() {}

bool Parametrization::operator==(const Parametrization& parametrization) const
{
    if (typeid(*this) != typeid(parametrization))
        return false;
    else
        return this->compare(parametrization);
}

bool Parametrization::operator!=(const Parametrization& parametrization) const
{
    return !(*this == parametrization);
}

bool Parametrization::compare(const Parametrization& parametrization) const
{
    if (particle_charge_ != parametrization.particle_charge_
        or particle_mass_ != parametrization.particle_mass_)
    if (lower_energy_lim_ != parametrization.lower_energy_lim_)
        return false;
    if (*medium_ != *parametrization.medium_)
        return false;
    else if (component_index_ != parametrization.component_index_)
        return false;
    else if (multiplier_ != parametrization.multiplier_)
        return false;
    else
        return true;
}

namespace PROPOSAL {
std::ostream& operator<<(std::ostream& os, Parametrization const& param)
{
    std::stringstream ss;
    ss << " Parametrization (" << &param << ") ";
    os << Helper::Centered(80, ss.str()) << '\n';

    os << "name: " << param.GetName() << '\n';

    param.print(os);

    os << "multiplier: " << param.multiplier_ << '\n';
    os << "current component index: " << param.component_index_ << '\n';
    os << "particle_mass: " << param.particle_mass_ << '\n';
    os << "particle_charge: " << param.particle_charge_ << '\n';
    os << "lower_energy_lim_: " << param.lower_energy_lim_ << '\n';
    os << *param.medium_ << '\n';
    os << Helper::Centered(80, "");
    return os;
}
} // namespace PROPOSAL

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

//----------------------------------------------------------------------------//
double Parametrization::FunctionToDNdxIntegral(double energy, double variable)
{
    return DifferentialCrossSection(energy, variable);
}

// ------------------------------------------------------------------------- //
double Parametrization::FunctionToDEdxIntegral(double energy, double variable)
{
    return variable * DifferentialCrossSection(energy, variable);
}

// ------------------------------------------------------------------------- //
double Parametrization::FunctionToDE2dxIntegral(double energy, double variable)
{
    return variable * variable * DifferentialCrossSection(energy, variable);
}

double Parametrization_builder::DifferentialCrossSection(
    double energy, double v)
{
    return diff_cross(*this, energy, v);
}

Parametrization::KinematicLimits Parametrization_builder::GetKinematicLimits(double energy)
{
    if (kinematic_limits)
        return kinematic_limits(*this, energy);

    Parametrization::KinematicLimits lim;
    lim.vMin = 0;
    lim.vMax = 0;

    return lim;
}

// ------------------------------------------------------------------------- //
// Getter
// ------------------------------------------------------------------------- //

size_t Parametrization::GetHash() const
{
    size_t seed = 0;
    hash_combine(seed, GetName(), std::abs(particle_charge_), particle_mass_,
        medium_->GetName());

    return seed;
}

size_t GetHash(const std::vector<Parametrization*>& params)
{
    size_t hash_digest = 0;
    for (auto param : params)
        hash_combine(hash_digest, param->GetHash());
    return hash_digest;
}
