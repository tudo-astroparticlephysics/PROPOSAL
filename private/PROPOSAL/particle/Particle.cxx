/*! \file   Particle.cxx
 *   \brief  Source file for the Particle routines.
 *
 *   For more details see the class documentation.
 *
 *   \date   2013.03.14
 *   \author Jan-Hendrik Koehne
 */

#include <cmath>
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

/******************************************************************************
 *                              Dynamic Particle                              *
 ******************************************************************************/

DynamicData::DynamicData(DynamicData::Type type)
    : type_id_(type)
    , position_(Vector3D())
    , direction_(Vector3D())
    , energy_(0)
    , parent_particle_energy_(0)
    , time_(0)
    , propagated_distance_(0)
{
}

DynamicData::DynamicData(const DynamicData& data)
    : type_id_(data.type_id_)
    , position_(data.position_)
    , direction_(data.direction_)
    , energy_(data.energy_)
    , parent_particle_energy_(data.parent_particle_energy_)
    , time_(data.time_)
    , propagated_distance_(data.propagated_distance_)
{
}

DynamicData::~DynamicData() {}

// ------------------------------------------------------------------------- //
std::ostream& PROPOSAL::operator<<(std::ostream& os, DynamicData const& data)
{
    std::stringstream ss;
    ss << " DynamicData (" << &data << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    os << "type: " << DynamicData::GetNameFromType(data.type_id_) << '\n';
    os << "position:" << '\n';
    os << data.position_ << '\n';
    os << "direction:" << '\n';
    os << data.direction_ << '\n';
    os << "energy: " << data.energy_ << '\n';
    os << "parent_particle_energy: " << data.parent_particle_energy_ << '\n';
    os << "time: " << data.time_ << '\n';
    os << "propagated distance: " << data.propagated_distance_ << '\n';

    data.print(os);

    os << Helper::Centered(60, "");
    return os;
}

// ------------------------------------------------------------------------- //
std::string DynamicData::GetNameFromType(Type type)
{
    switch (type)
    {
        case None:
            return "None";
        case Particle:
            return "Particle";
        case Brems:
            return "Brems";
        case DeltaE:
            return "DeltaE";
        case Epair:
            return "Epair";
        case NuclInt:
            return "NuclInt";
        case MuPair:
            return "MuPair";
        case Hadrons:
            return "Hadrons";
        case ContinuousEnergyLoss:
            return "ContinuousEnergyLoss";
        case WeakInt:
            return "WeakInt";
        default:
            return "Type not found";
    }
}

/******************************************************************************
 *                              Particle                                       *
 ******************************************************************************/

Particle::Particle()
    : DynamicData(DynamicData::Particle)
    , particle_def_(MuMinusDef::Get())
    , momentum_(0)
    , parent_particle_id_(0)
    , particle_id_(1)
    , entry_point_(Vector3D())
    , entry_time_(0)
    , entry_energy_(0)
    , exit_point_(Vector3D())
    , exit_time_(0)
    , exit_energy_(0)
    , closest_approach_point_(Vector3D())
    , closest_approach_time_(0)
    , closest_approach_energy_(0)
    , elost_(0)
{
    SetEnergy(energy_);
}

Particle::Particle(const Particle& particle)
    : DynamicData(particle)
    , particle_def_(particle.particle_def_)
    , momentum_(particle.momentum_)
    , parent_particle_id_(particle.parent_particle_id_)
    , particle_id_(particle.particle_id_)
    , entry_point_(particle.entry_point_)
    , entry_time_(particle.entry_time_)
    , entry_energy_(particle.entry_energy_)
    , exit_point_(particle.exit_point_)
    , exit_time_(particle.exit_time_)
    , exit_energy_(particle.exit_energy_)
    , closest_approach_point_(particle.closest_approach_point_)
    , closest_approach_time_(particle.closest_approach_time_)
    , closest_approach_energy_(particle.closest_approach_energy_)
    , elost_(particle.elost_)
{
}

Particle::Particle(const ParticleDef& particleDef)
    : DynamicData(DynamicData::Particle)
    , particle_def_(particleDef)
    , momentum_(0)
    , parent_particle_id_(0)
    , particle_id_(1)
    , entry_point_(Vector3D())
    , entry_time_(0)
    , entry_energy_(0)
    , exit_point_(Vector3D())
    , exit_time_(0)
    , exit_energy_(0)
    , closest_approach_point_(Vector3D())
    , closest_approach_time_(0)
    , closest_approach_energy_(0)
    , elost_(0)
{
    SetEnergy(energy_);
}

// ------------------------------------------------------------------------- //
// Operators
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
bool Particle::operator==(const Particle& particle) const
{
    if (propagated_distance_ != particle.propagated_distance_)
        return false;
    if (position_ != particle.position_)
        return false;
    if (direction_ != particle.direction_)
        return false;
    if (momentum_ != particle.momentum_)
        return false;
    if (energy_ != particle.energy_)
        return false;
    if (parent_particle_id_ != particle.parent_particle_id_)
        return false;
    if (parent_particle_energy_ != particle.parent_particle_energy_)
        return false;
    if (particle_id_ != particle.particle_id_)
        return false;
    if (entry_point_ != particle.entry_point_)
        return false;
    if (entry_time_ != particle.entry_time_)
        return false;
    if (entry_energy_ != particle.entry_energy_)
        return false;
    if (exit_point_ != particle.exit_point_)
        return false;
    if (exit_time_ != particle.exit_time_)
        return false;
    if (exit_energy_ != particle.exit_energy_)
        return false;
    if (closest_approach_point_ != particle.closest_approach_point_)
        return false;
    if (closest_approach_time_ != particle.closest_approach_time_)
        return false;
    if (closest_approach_energy_ != particle.closest_approach_energy_)
        return false;
    if (elost_ != particle.elost_)
        return false;
    if (particle_def_ != particle.particle_def_)
    {
        return false;
    }

    return true;
}

// ------------------------------------------------------------------------- //
bool Particle::operator!=(const Particle& particle) const
{
    return !(*this == particle);
}

// ------------------------------------------------------------------------- //
// Methods
// ------------------------------------------------------------------------- //

void Particle::InjectState(const Particle& particle)
{
    position_                = particle.position_;
    direction_               = particle.direction_;
    energy_                  = particle.energy_;
    parent_particle_energy_  = particle.parent_particle_energy_;
    time_                    = particle.time_;
    propagated_distance_     = particle.propagated_distance_;
    momentum_                = particle.momentum_;
    entry_point_             = particle.entry_point_;
    entry_time_              = particle.entry_time_;
    entry_energy_            = particle.entry_energy_;
    exit_point_              = particle.exit_point_;
    exit_time_               = particle.exit_time_;
    exit_energy_             = particle.exit_energy_;
    closest_approach_point_  = particle.closest_approach_point_;
    closest_approach_time_   = particle.closest_approach_time_;
    closest_approach_energy_ = particle.closest_approach_energy_;
    elost_                   = particle.elost_;
}

// ------------------------------------------------------------------------- //
// Setter
// ------------------------------------------------------------------------- //

void Particle::SetEnergy(double energy)
{
    energy_   = std::max(energy, particle_def_.mass);
    momentum_ = std::sqrt(std::max((energy_ + particle_def_.mass) * (energy_ - particle_def_.mass), 0.0));
}

// ------------------------------------------------------------------------- //
void Particle::SetMomentum(double momentum)
{
    momentum_ = momentum;
    energy_   = std::sqrt(momentum_ * momentum_ + particle_def_.mass * particle_def_.mass);
}

// ------------------------------------------------------------------------- //
// Print
// ------------------------------------------------------------------------- //

void Particle::print(std::ostream& os) const
{
    os << "definition:" << '\n';
    os << particle_def_ << '\n';
    os << "momentum [MeV]: " << momentum_ << '\n';
    os << "energy lost in detector [MeV]: " << elost_ << '\n';

    os << "detector entry point:" << '\n';
    os << entry_point_ << '\n';
    os << "entry time [s]: " << entry_time_ << '\n';
    os << "entry energy [MeV]: " << entry_energy_ << '\n';
    os << "detector exit point:" << '\n';
    os << exit_point_ << '\n';
    os << "exit time [s]: " << exit_time_ << '\n';
    os << "exit energy [MeV]: " << exit_energy_ << '\n';
    os << "detector closest approach point:" << '\n';
    os << closest_approach_point_ << '\n';
    os << "closest approach time [s]: " << closest_approach_time_ << '\n';
    os << "closest approach energy [MeV]: " << closest_approach_energy_ << '\n';
}
