#include "PROPOSAL/secondaries/compton/NaivCompton.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>
#include <stdexcept>

using std::fmod;
using std::get;
using std::make_tuple;
using std::sqrt;

using namespace PROPOSAL;
double secondaries::NaivCompton::CalculateRho(
    double primary_energy, double loss_energy)
{
    return loss_energy / primary_energy;
}

tuple<Vector3D, Vector3D> secondaries::NaivCompton::CalculateDirections(
    Vector3D primary_dir, double energy, double v, double rnd)
{
    auto com_energy = energy + ME; // center of mass energy
    auto kin_energy = energy;      // kinetic energy
    auto cosphi_gamma = (com_energy * (1. - v) - ME)
        / ((1. - v) * sqrt(com_energy * kin_energy));
    auto cosphi_electron
        = (com_energy * v - ME) / (v * sqrt(com_energy * kin_energy));
    auto rnd_theta = rnd * 2. * PI;
    auto dir_gamma = deflect(primary_dir, cosphi_gamma, rnd_theta);
    auto dir_electron
        = deflect(primary_dir, cosphi_electron, fmod(rnd_theta + PI, 2. * PI));
    return make_tuple(dir_gamma, dir_electron);
}

tuple<double, double> secondaries::NaivCompton::CalculateEnergy(
    double energy, double v)
{
    return make_tuple(energy * (1 - v), energy * v);
}

vector<ParticleState> secondaries::NaivCompton::CalculateSecondaries(
        StochasticLoss loss, const Component&, vector<double> &rnd)
{
    auto v = loss.loss_energy /  loss.parent_particle_energy;
    auto secondary_energies = CalculateEnergy(loss.parent_particle_energy, v);
    auto secondary_dir = CalculateDirections(loss.direction, loss.loss_energy,
                                             v, rnd[1]);

    auto sec = std::vector<ParticleState>();
    sec.emplace_back(ParticleType::Gamma, loss.position,
                     get<GAMMA>(secondary_dir), get<GAMMA>(secondary_energies),
                             loss.time, 0.);
    sec.emplace_back(ParticleType::EMinus, loss.position,
                     get<EMINUS>(secondary_dir), get<EMINUS>(secondary_energies),
                             loss.time, 0.);
    return sec;
}
