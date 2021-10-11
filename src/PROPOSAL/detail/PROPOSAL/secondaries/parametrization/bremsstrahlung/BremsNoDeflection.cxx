#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/BremsNoDeflection.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>
#include <stdexcept>

using namespace PROPOSAL;

std::vector<ParticleState> secondaries::BremsNoDeflection::CalculateSecondaries(
    StochasticLoss loss, const Component& c, std::vector<double>& rnd)
{
    auto directions = CalculateDirections(
            loss.direction, loss.parent_particle_energy, loss.energy, c, rnd);

    auto primary_lepton = ParticleState();
    primary_lepton.energy = loss.parent_particle_energy - loss.energy;
    primary_lepton.type = primary_lepton_type;
    primary_lepton.time = loss.time;
    primary_lepton.position = loss.position;
    primary_lepton.direction = directions.first;
    primary_lepton.propagated_distance = 0.;

    auto brems_photon = ParticleState();
    brems_photon.energy = loss.energy;
    brems_photon.SetType(PROPOSAL::ParticleType::Gamma);
    brems_photon.time = loss.time;
    brems_photon.position = loss.position;
    brems_photon.direction = directions.second;

    auto sec = std::vector<ParticleState>{primary_lepton, brems_photon};
    return sec;
}

std::pair<Cartesian3D, Cartesian3D>
        secondaries::BremsNoDeflection::CalculateDirections(
                const Vector3D& init_direction, double, double,
                const Component&, std::vector<double>&) {
    return std::pair<Cartesian3D, Cartesian3D>(init_direction, init_direction);
}