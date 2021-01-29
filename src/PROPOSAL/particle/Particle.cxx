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
#include <sstream>
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
    , position(Cartesian3D())
    , direction(Cartesian3D())
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

StochasticLoss::StochasticLoss(int type, double loss_energy, const Vector3D& position,
                               const Vector3D& direction, double time,
                               double propagated_distance,
                               double parent_particle_energy)
                               : Loss(type, loss_energy, parent_particle_energy),
                               position(position), direction(direction), time(time),
                               propagated_distance(propagated_distance) {}

ContinuousLoss::ContinuousLoss(double energy, double parent_particle_energy,
                               const Vector3D& start_position, double length,
                               const Vector3D& direction_initial,
                               const Vector3D& direction_final,
                               double time_initial, double time_final)
                               : Loss((int)InteractionType::ContinuousEnergyLoss,
                                      energy, parent_particle_energy),
                                      start_position(start_position),
                                      length(length),
                                      direction_initial(direction_initial),
                                      direction_final(direction_final),
                                      time_initial(time_initial),
                                      time_final(time_final) {}