#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/NaivBremsstrahlung.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>
#include <stdexcept>

using namespace PROPOSAL;

std::vector<ParticleState> secondaries::NaivBremsstrahlung::CalculateSecondaries(
    StochasticLoss loss, const Component&, std::vector<double>&)
{
    auto primary_lepton = ParticleState();
    primary_lepton.energy = loss.parent_particle_energy - loss.energy;
    primary_lepton.type = primary_lepton_type;
    primary_lepton.time = loss.time;
    primary_lepton.position = loss.position;
    primary_lepton.direction = loss.direction;
    primary_lepton.propagated_distance = 0.;

    auto brems_photon = ParticleState();
    brems_photon.energy = loss.energy;
    brems_photon.SetType(PROPOSAL::ParticleType::Gamma);
    brems_photon.time = loss.time;
    brems_photon.position = loss.position;
    brems_photon.direction = loss.direction;

    auto sec = std::vector<ParticleState>{primary_lepton, brems_photon};
    return sec;
}
