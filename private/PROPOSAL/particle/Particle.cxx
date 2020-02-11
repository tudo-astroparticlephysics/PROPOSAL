/*! \file   Particle.cxx
 *   \brief  Source file for the Particle routines.
 *
 *   For more details see the class documentation.
 *
 *   \date   2013.03.14
 *   \author Jan-Hendrik Koehne
 */

#include <cmath>
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/methods.h"
#include <iostream>
#include <string>

using namespace PROPOSAL;


/* std::map<const int, const ParticleDef&> Id_Particle_Map; */
/* std::map<const int, std::string> Id_Interaction_Name_Map; */

/******************************************************************************
 *                              Dynamic Particle                              *
 ******************************************************************************/

DynamicData::DynamicData()
    : type_id_(0)
    , position_(Vector3D())
    , direction_(Vector3D())
    , energy_(0)
    , parent_particle_energy_(0)
    , time_(0)
    , propagated_distance_(0)
{
}

DynamicData::DynamicData(const int& type)
    : type_id_(type)
    , position_(Vector3D())
    , direction_(Vector3D())
    , energy_(0)
    , parent_particle_energy_(0)
    , time_(0)
    , propagated_distance_(0)
{
}

DynamicData::DynamicData(const int& type, const Vector3D& position, const Vector3D& direction, const double& energy, const double& parent_particle_energy, const double& time, const double& distance)
    : type_id_(type)
    , position_(position)
    , direction_(direction)
    , energy_(energy)
    , parent_particle_energy_(parent_particle_energy)
    , time_(time)
    , propagated_distance_(distance)
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



DynamicData::DynamicData(DynamicData&& other)
    : type_id_(std::move(other.type_id_))
    , position_(std::move(other.position_))
    , direction_(std::move(other.direction_))
    , energy_(std::move(other.energy_))
    , parent_particle_energy_(std::move(other.parent_particle_energy_))
    , time_(std::move(other.time_))
    , propagated_distance_(std::move(other.propagated_distance_))
{
}

DynamicData::~DynamicData() {}

DynamicData& DynamicData::operator=(const DynamicData& data)
{
    if (this != &data){
        type_id_ = data.type_id_;
        position_ = data.position_;
        direction_ = data.direction_;
        energy_ = data.energy_;
        parent_particle_energy_ = data.parent_particle_energy_;
        time_ = data.time_;
        propagated_distance_ = data.propagated_distance_;
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool DynamicData::operator==(const DynamicData& dynamic_data) const
{
    if (type_id_ != dynamic_data.type_id_)
        return false;
    if (position_ != dynamic_data.position_)
        return false;
    if (direction_ != dynamic_data.direction_)
        return false;
    if (energy_ != dynamic_data.energy_)
        return false;
    if (parent_particle_energy_ != dynamic_data.parent_particle_energy_)
        return false;
    if (time_ != dynamic_data.time_)
        return false;
    if (propagated_distance_ != dynamic_data.propagated_distance_)
        return false;

    return true;
}

// ------------------------------------------------------------------------- //
bool DynamicData::operator!=(const DynamicData& dynamic_data) const
{
    return !(*this == dynamic_data);
}

// ------------------------------------------------------------------------- //
std::ostream& PROPOSAL::operator<<(std::ostream& os, DynamicData const& data)
{
    std::stringstream ss;
    ss << " DynamicData (" << &data << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    os << "type: " << data.GetName() << '\n';
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
std::string DynamicData::GetName() const
{
    auto p_search = Id_Particle_Map.find(type_id_);
    if (p_search != Id_Particle_Map.end()) {
        return p_search->second.name;
    }

    auto i_search = Id_Interaction_Name_Map.find(type_id_);
    if (i_search != Id_Interaction_Name_Map.end()) {
        return i_search->second;
    }

    return "Not found." ;
}

void DynamicData::SetEnergy(double energy)
{
    auto particle = Id_Particle_Map.find(type_id_);

    if (particle != Id_Particle_Map.end()) {
        energy_   = std::max(energy, particle->second.mass);
    } else {
        energy_ = energy;
    }
}

// ------------------------------------------------------------------------- //
void DynamicData::SetMomentum(double momentum)
{
    auto particle = Id_Particle_Map.find(type_id_);

    if (particle != Id_Particle_Map.end()) {
        energy_   = std::sqrt(momentum * momentum + particle->second.mass * particle->second.mass);
    } else {
        energy_ = momentum;
    }
}



double DynamicData::GetMomentum() const
{
    auto particle = Id_Particle_Map.find(type_id_);

    if (particle != Id_Particle_Map.end()) {
        return std::sqrt((energy_ + particle->second.mass) * (energy_ - particle->second.mass));
    }

    return energy_;

}

void DynamicData::DeflectDirection(double cosphi_deflect, double theta_deflect) {

    Vector3D old_direction = GetDirection();

    old_direction.CalculateSphericalCoordinates();
    double sinphi_deflect = std::sqrt( std::max(0., (1. - cosphi_deflect) * (1. + cosphi_deflect) ));
    double tx = sinphi_deflect * std::cos(theta_deflect);
    double ty = sinphi_deflect * std::sin(theta_deflect);
    double tz = std::sqrt(std::max(1. - tx * tx - ty * ty, 0.));
    if(cosphi_deflect < 0. ){
        // Backward deflection
        tz = -tz;
    }

    long double sinth, costh, sinph, cosph;
    sinth = (long double)std::sin(old_direction.GetTheta());
    costh = (long double)std::cos(old_direction.GetTheta());
    sinph = (long double)std::sin(old_direction.GetPhi());
    cosph = (long double)std::cos(old_direction.GetPhi());

    const Vector3D rotate_vector_x = Vector3D(costh * cosph, costh * sinph, -sinth);
    const Vector3D rotate_vector_y = Vector3D(-sinph, cosph, 0.);

    // Rotation towards all tree axes
    Vector3D new_direction( tz * old_direction + tx * rotate_vector_x + ty * rotate_vector_y );

    direction_ = new_direction;
    direction_.CalculateSphericalCoordinates();
}

/******************************************************************************
 *                              Particle                                       *
 ******************************************************************************/

// Particle::Particle()
//     : DynamicData(static_cast<int>(InteractionType::Particle))
//     , particle_def_(MuMinusDef::Get())
//     , momentum_(0)
//     , parent_particle_id_(0)
//     , particle_id_(1)
//     , entry_point_(Vector3D())
//     , entry_time_(0)
//     , entry_energy_(0)
//     , exit_point_(Vector3D())
//     , exit_time_(0)
//     , exit_energy_(0)
//     , closest_approach_point_(Vector3D())
//     , closest_approach_time_(0)
//     , closest_approach_energy_(0)
//     , elost_(0)
// {
//     SetEnergy(energy_);
// }

// Particle::Particle(const Particle& particle)
//     : DynamicData(particle)
//     , particle_def_(particle.particle_def_)
//     , momentum_(particle.momentum_)
//     , parent_particle_id_(particle.parent_particle_id_)
//     , particle_id_(particle.particle_id_)
//     , entry_point_(particle.entry_point_)
//     , entry_time_(particle.entry_time_)
//     , entry_energy_(particle.entry_energy_)
//     , exit_point_(particle.exit_point_)
//     , exit_time_(particle.exit_time_)
//     , exit_energy_(particle.exit_energy_)
//     , closest_approach_point_(particle.closest_approach_point_)
//     , closest_approach_time_(particle.closest_approach_time_)
//     , closest_approach_energy_(particle.closest_approach_energy_)
//     , elost_(particle.elost_)
// {
// }

// Particle::Particle(const ParticleDef& particleDef)
//     : DynamicData(static_cast<int>(InteractionType::Particle))
//     , particle_def_(particleDef)
//     , momentum_(0)
//     , parent_particle_id_(0)
//     , particle_id_(1)
//     , entry_point_(Vector3D())
//     , entry_time_(0)
//     , entry_energy_(0)
//     , exit_point_(Vector3D())
//     , exit_time_(0)
//     , exit_energy_(0)
//     , closest_approach_point_(Vector3D())
//     , closest_approach_time_(0)
//     , closest_approach_energy_(0)
//     , elost_(0)
// {
//     SetEnergy(energy_);
// }

// // ------------------------------------------------------------------------- //
// // Operators
// // ------------------------------------------------------------------------- //

// // ------------------------------------------------------------------------- //
// bool Particle::operator==(const Particle& particle) const
// {
//     if (propagated_distance_ != particle.propagated_distance_)
//         return false;
//     if (position_ != particle.position_)
//         return false;
//     if (direction_ != particle.direction_)
//         return false;
//     if (momentum_ != particle.momentum_)
//         return false;
//     if (energy_ != particle.energy_)
//         return false;
//     if (parent_particle_id_ != particle.parent_particle_id_)
//         return false;
//     if (parent_particle_energy_ != particle.parent_particle_energy_)
//         return false;
//     if (particle_id_ != particle.particle_id_)
//         return false;
//     if (entry_point_ != particle.entry_point_)
//         return false;
//     if (entry_time_ != particle.entry_time_)
//         return false;
//     if (entry_energy_ != particle.entry_energy_)
//         return false;
//     if (exit_point_ != particle.exit_point_)
//         return false;
//     if (exit_time_ != particle.exit_time_)
//         return false;
//     if (exit_energy_ != particle.exit_energy_)
//         return false;
//     if (closest_approach_point_ != particle.closest_approach_point_)
//         return false;
//     if (closest_approach_time_ != particle.closest_approach_time_)
//         return false;
//     if (closest_approach_energy_ != particle.closest_approach_energy_)
//         return false;
//     if (elost_ != particle.elost_)
//         return false;
//     if (particle_def_ != particle.particle_def_)
//     {
//         return false;
//     }

//     return true;
// }

// // ------------------------------------------------------------------------- //
// bool Particle::operator!=(const Particle& particle) const
// {
//     return !(*this == particle);
// }

// Particle Particle::operator=(const Particle&)
// {
//     return Particle(*this);
// }

// // ------------------------------------------------------------------------- //
// // Methods
// // ------------------------------------------------------------------------- //

// void Particle::InjectState(const Particle& particle)
// {
//     position_                = particle.position_;
//     direction_               = particle.direction_;
//     energy_                  = particle.energy_;
//     parent_particle_energy_  = particle.parent_particle_energy_;
//     time_                    = particle.time_;
//     propagated_distance_     = particle.propagated_distance_;
//     momentum_                = particle.momentum_;
//     entry_point_             = particle.entry_point_;
//     entry_time_              = particle.entry_time_;
//     entry_energy_            = particle.entry_energy_;
//     exit_point_              = particle.exit_point_;
//     exit_time_               = particle.exit_time_;
//     exit_energy_             = particle.exit_energy_;
//     closest_approach_point_  = particle.closest_approach_point_;
//     closest_approach_time_   = particle.closest_approach_time_;
//     closest_approach_energy_ = particle.closest_approach_energy_;
//     elost_                   = particle.elost_;
// }

// // ------------------------------------------------------------------------- //
// // Setter
// // ------------------------------------------------------------------------- //

// /* void Particle::SetEnergy(double energy) */
// /* { */
// /*     energy_   = std::max(energy, particle_def_.mass); */
// /*     momentum_ = std::sqrt(std::max((energy_ + particle_def_.mass) * (energy_ - particle_def_.mass), 0.0)); */
// /* } */

// // ------------------------------------------------------------------------- //
// /* void Particle::SetMomentum(double momentum) */
// /* { */
// /*     momentum_ = momentum; */
// /*     energy_   = std::sqrt(momentum_ * momentum_ + particle_def_.mass * particle_def_.mass); */
// /* } */

// void Particle::DeflectDirection(double cosphi_deflect, double theta_deflect) {

//     Vector3D old_direction = GetDirection();

//     old_direction.CalculateSphericalCoordinates();
//     double sinphi_deflect = std::sqrt( std::max(0., (1. - cosphi_deflect) * (1. + cosphi_deflect) ));
//     double tx = sinphi_deflect * std::cos(theta_deflect);
//     double ty = sinphi_deflect * std::sin(theta_deflect);
//     double tz = std::sqrt(std::max(1. - tx * tx - ty * ty, 0.));
//     if(cosphi_deflect < 0. ){
//         // Backward deflection
//         tz = -tz;
//     }

//     long double sinth, costh, sinph, cosph;
//     sinth = (long double)std::sin(old_direction.GetTheta());
//     costh = (long double)std::cos(old_direction.GetTheta());
//     sinph = (long double)std::sin(old_direction.GetPhi());
//     cosph = (long double)std::cos(old_direction.GetPhi());

//     const Vector3D rotate_vector_x = Vector3D(costh * cosph, costh * sinph, -sinth);
//     const Vector3D rotate_vector_y = Vector3D(-sinph, cosph, 0.);

//     // Rotation towards all tree axes
//     Vector3D new_direction( tz * old_direction + tx * rotate_vector_x + ty * rotate_vector_y );

//     SetDirection(new_direction);

// }

// // ------------------------------------------------------------------------- //
// // Print
// // ------------------------------------------------------------------------- //

// void Particle::print(std::ostream& os) const
// {
//     os << "definition:" << '\n';
//     os << particle_def_ << '\n';
//     os << "momentum [MeV]: " << momentum_ << '\n';
//     os << "energy lost in detector [MeV]: " << elost_ << '\n';

//     os << "detector entry point:" << '\n';
//     os << entry_point_ << '\n';
//     os << "entry time [s]: " << entry_time_ << '\n';
//     os << "entry energy [MeV]: " << entry_energy_ << '\n';
//     os << "detector exit point:" << '\n';
//     os << exit_point_ << '\n';
//     os << "exit time [s]: " << exit_time_ << '\n';
//     os << "exit energy [MeV]: " << exit_energy_ << '\n';
//     os << "detector closest approach point:" << '\n';
//     os << closest_approach_point_ << '\n';
//     os << "closest approach time [s]: " << closest_approach_time_ << '\n';
//     os << "closest approach energy [MeV]: " << closest_approach_energy_ << '\n';
// }
