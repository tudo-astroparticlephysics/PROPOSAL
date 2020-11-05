/*! \file   Particle.cxx
 *   \brief  Source file for the Particle routines.
 *
 *   For more details see the class documentation.
 *
 *   \date   2013.03.14
 *   \author Jan-Hendrik Koehne
 */

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include <cmath>
#include <iostream>
#include <string>

using namespace PROPOSAL;

namespace PROPOSAL {

std::ostream& operator<<(std::ostream& os, ParticleState const& data)
{
    std::stringstream ss;
    ss << " ParticleState (" << &data << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    os << "type: " << data.type << '\n';
    os << "position:" << '\n';
    os << data.position << '\n';
    os << "direction:" << '\n';
    os << data.direction << '\n';
    os << "energy: " << data.energy << '\n';
    os << "time: " << data.time << '\n';
    os << "propagated distance: " << data.propagated_distance << '\n';

    data.print(os);

    os << Helper::Centered(60, "");
    return os;
}

} // namespace PROPOSAL

ParticleState::ParticleState()
    : type(0)
    , position(Vector3D())
    , direction(Vector3D())
    , energy(0)
    , time(0)
    , propagated_distance(0)
{
}

ParticleState::ParticleState(const Vector3D& position, const Vector3D& direction,
                             const double& energy, const double& time,
                             const double& distance)
    : type(static_cast<int>(ParticleType::None))
    , position(position)
    , direction(direction)
    , energy(energy)
    , time(time)
    , propagated_distance(distance)
{
}

ParticleState::ParticleState(const ParticleType& type, const Vector3D& position,
                             const Vector3D& direction, const double& energy, const double& time,
                             const double& distance)
    : type(static_cast<int>(type))
    , position(position)
    , direction(direction)
    , energy(energy)
    , time(time)
    , propagated_distance(distance)
{
}

bool ParticleState::operator==(const ParticleState& dynamic_data) const
{
    if (type != dynamic_data.type)
        return false;
    if (position != dynamic_data.position)
        return false;
    if (direction != dynamic_data.direction)
        return false;
    if (energy != dynamic_data.energy)
        return false;
    if (time != dynamic_data.time)
        return false;
    if (propagated_distance != dynamic_data.propagated_distance)
        return false;

    return true;
}

bool ParticleState::operator!=(const ParticleState& dynamic_data) const
{
    return !(*this == dynamic_data);
}

ParticleDef ParticleState::GetParticleDef() const
{
    auto p_search = Type_Particle_Map.find(static_cast<ParticleType>(type));
    if (p_search != Type_Particle_Map.end()) {
        return p_search->second;
    }

    throw std::invalid_argument("Particle def for ParticleType not found.");
}

void ParticleState::SetMomentum(double momentum)
{
    auto particle = Type_Particle_Map.find(static_cast<ParticleType>(type));
    if (particle != Type_Particle_Map.end())
        energy = std::sqrt(momentum * momentum
            + particle->second.mass * particle->second.mass);
    else
        energy = momentum;
}

double ParticleState::GetMomentum() const
{
    auto particle = Type_Particle_Map.find(static_cast<ParticleType>(type));
    if (particle != Type_Particle_Map.end())
        return std::sqrt((energy + particle->second.mass)
            * (energy - particle->second.mass));
    return energy;
}

void ParticleState::DeflectDirection(double cosphi_deflect, double theta_deflect)
{

    auto old_direction = direction;

    old_direction.CalculateSphericalCoordinates();
    auto sinphi_deflect = std::sqrt(
        std::max(0., (1. - cosphi_deflect) * (1. + cosphi_deflect)));
    auto tx = sinphi_deflect * std::cos(theta_deflect);
    auto ty = sinphi_deflect * std::sin(theta_deflect);
    auto tz = std::sqrt(std::max(1. - tx * tx - ty * ty, 0.));
    if (cosphi_deflect < 0.)
        tz = -tz; // Backward deflection

    auto sinth = std::sin(old_direction.GetTheta());
    auto costh = std::cos(old_direction.GetTheta());
    auto sinph = std::sin(old_direction.GetPhi());
    auto cosph = std::cos(old_direction.GetPhi());

    auto rotate_vector_x = Vector3D(costh * cosph, costh * sinph, -sinth);
    auto rotate_vector_y = Vector3D(-sinph, cosph, 0.);

    // Rotation towards all tree axes
    auto new_direction = Vector3D(
        tz * old_direction + tx * rotate_vector_x + ty * rotate_vector_y);
    new_direction.CalculateSphericalCoordinates();

    direction = new_direction;
}

StochasticLoss::StochasticLoss(int type, double loss_energy, Vector3D position,
                               Vector3D direction, double time,
                               double propagated_distance,
                               double parent_particle_energy) : Loss(type, loss_energy),
                               position(position), direction(direction), time(time),
                               propagated_distance(propagated_distance),
                               parent_particle_energy(parent_particle_energy) {}

ContinuousLoss::ContinuousLoss(std::pair<double, double> energies,
                               std::pair<Vector3D, Vector3D> positions,
                               std::pair<Vector3D, Vector3D> directions,
                               std::pair<double, double> times)
                               : Loss((int)InteractionType::ContinuousEnergyLoss,
                                      std::abs(energies.second - energies.first)),
                                      energies(energies),
                                      positions(positions),
                                      directions(directions),
                                      times(times) {}