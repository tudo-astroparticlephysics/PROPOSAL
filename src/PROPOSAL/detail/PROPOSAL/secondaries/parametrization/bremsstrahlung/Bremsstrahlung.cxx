#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/Bremsstrahlung.h"
#include "PROPOSAL/particle/Particle.h"

using namespace PROPOSAL;

std::vector<ParticleState> secondaries::Bremsstrahlung::CalculateSecondaries(
        StochasticLoss loss, const Component& c, std::vector<double>& rnd)
        {
    auto directions = CalculateDirections(
            loss.direction, loss.parent_particle_energy, loss.energy, c, rnd);

    ParticleState primary_lepton;
    primary_lepton.energy = loss.parent_particle_energy - loss.energy;
    primary_lepton.type = primary_lepton_type;
    primary_lepton.time = loss.time;
    primary_lepton.position = loss.position;
    primary_lepton.direction = directions.first;
    primary_lepton.propagated_distance = 0.;

    ParticleState brems_photon;
    brems_photon.energy = loss.energy;
    brems_photon.SetType(PROPOSAL::ParticleType::Gamma);
    brems_photon.time = loss.time;
    brems_photon.position = loss.position;
    brems_photon.direction = directions.second;

    auto sec = std::vector<ParticleState>{primary_lepton, brems_photon};
    return sec;
        }