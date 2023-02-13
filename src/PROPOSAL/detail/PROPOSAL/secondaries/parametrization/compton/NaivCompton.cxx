#include "PROPOSAL/secondaries/parametrization/compton/NaivCompton.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>
#include <stdexcept>

using namespace PROPOSAL;
double secondaries::NaivCompton::CalculateRho(
    double primary_energy, double loss_energy)
{
    return loss_energy / primary_energy;
}

std::tuple<Cartesian3D, Cartesian3D> secondaries::NaivCompton::CalculateDirections(
    const Vector3D& primary_dir, double energy, double v, double rnd)
{
    // use energy-momentum conservation formula
    auto cosphi_gamma = 1 - (v * ME) / (energy * (1 - v));
    auto cosphi_electron = v * (energy + ME) / std::sqrt(2 * v * energy * ME + v * v * energy * energy);

    auto rnd_theta = rnd * 2. * PI;
    auto dir_gamma = Cartesian3D(primary_dir);
    dir_gamma.deflect(cosphi_gamma, rnd_theta);
    auto dir_electron = Cartesian3D(primary_dir);
    dir_electron.deflect(cosphi_electron, std::fmod(rnd_theta + PI, 2. * PI));
    return std::make_tuple(dir_gamma, dir_electron);
}

std::tuple<double, double> secondaries::NaivCompton::CalculateEnergy(
    double energy, double v)
{
    return std::make_tuple(energy * (1 - v), ME + energy * v);
}

std::vector<ParticleState> secondaries::NaivCompton::CalculateSecondaries(
        StochasticLoss loss, const Component&, std::vector<double> &rnd)
{
    auto v = loss.energy /  loss.parent_particle_energy;
    auto secondary_energies = CalculateEnergy(loss.parent_particle_energy, v);
    auto secondary_dir = CalculateDirections(loss.direction, loss.parent_particle_energy, v,
                                             rnd[0]);

    auto sec = std::vector<ParticleState>();
    sec.emplace_back(ParticleType::Gamma, loss.position,
                     std::get<GAMMA>(secondary_dir), std::get<GAMMA>(secondary_energies),
                             loss.time, 0.);
    sec.emplace_back(ParticleType::EMinus, loss.position,
                     std::get<EMINUS>(secondary_dir), std::get<EMINUS>(secondary_energies),
                             loss.time, 0.);
    return sec;
}
