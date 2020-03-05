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

namespace PROPOSAL {

// ------------------------------------------------------------------------- //
std::ostream& operator<<(std::ostream& os, DynamicData const& data)
{
    std::stringstream ss;
    ss << " DynamicData (" << &data << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    os << "type: " << data.GetType() << '\n';
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

} // namespace PROPOSAL

/******************************************************************************
 *                              Dynamic Particle                              *
 ******************************************************************************/

DynamicData::DynamicData()
    : type_(0)
    , position_(Vector3D())
    , direction_(Vector3D())
    , energy_(0)
    , parent_particle_energy_(0)
    , time_(0)
    , propagated_distance_(0)
{
}

DynamicData::DynamicData(const int& type)
    : type_(type)
    , position_(Vector3D())
    , direction_(Vector3D())
    , energy_(0)
    , parent_particle_energy_(0)
    , time_(0)
    , propagated_distance_(0)
{
}

DynamicData::DynamicData(const int& type, const Vector3D& position, const Vector3D& direction, const double& energy, const double& parent_particle_energy, const double& time, const double& distance)
    : type_(type)
    , position_(position)
    , direction_(direction)
    , energy_(energy)
    , parent_particle_energy_(parent_particle_energy)
    , time_(time)
    , propagated_distance_(distance)
{
}

DynamicData::DynamicData(const DynamicData& data)
    : type_(data.type_)
    , position_(data.position_)
    , direction_(data.direction_)
    , energy_(data.energy_)
    , parent_particle_energy_(data.parent_particle_energy_)
    , time_(data.time_)
    , propagated_distance_(data.propagated_distance_)
{
}



DynamicData::DynamicData(DynamicData&& other)
    : type_(std::move(other.type_))
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
        type_ = data.type_;
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
    if (type_ != dynamic_data.type_)
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
std::string DynamicData::GetName() const
{
    auto p_search = Type_Particle_Map.find(type_);
    if (p_search != Type_Particle_Map.end()) {
        return p_search->second.name;
    }

    auto i_search = Type_Interaction_Name_Map.find(type_);
    if (i_search != Type_Interaction_Name_Map.end()) {
        return i_search->second;
    }

    return "Not found." ;
}

void DynamicData::SetEnergy(double energy)
{
    auto particle = Type_Particle_Map.find(type_);

    if (particle != Type_Particle_Map.end()) {
        energy_   = std::max(energy, particle->second.mass);
    } else {
        energy_ = energy;
    }
}

// ------------------------------------------------------------------------- //
void DynamicData::SetMomentum(double momentum)
{
    auto particle = Type_Particle_Map.find(type_);

    if (particle != Type_Particle_Map.end()) {
        energy_   = std::sqrt(momentum * momentum + particle->second.mass * particle->second.mass);
    } else {
        energy_ = momentum;
    }
}



double DynamicData::GetMomentum() const
{
    auto particle = Type_Particle_Map.find(type_);

    if (particle != Type_Particle_Map.end()) {
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
