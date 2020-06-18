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
};

tuple<Vector3D, Vector3D> secondaries::NaivCompton::CalculateDirections(
    Vector3D primary_dir, double energy, double rho, double rnd)
{
    auto com_energy = energy + ME; // center of mass energy
    auto kin_energy = energy;      // kinetic energy
    auto cosphi_gamma = (com_energy * (1. - rho) - ME)
        / ((1. - rho) * sqrt(com_energy * kin_energy));
    auto cosphi_electron
        = (com_energy * rho - ME) / (rho * sqrt(com_energy * kin_energy));
    auto rnd_theta = rnd * 2. * PI;
    auto dir_gamma = deflect(primary_dir, cosphi_gamma, rnd_theta);
    auto dir_electron
        = deflect(primary_dir, cosphi_electron, fmod(rnd_theta + PI, 2. * PI));
    return make_tuple(dir_gamma, dir_electron);
}

tuple<double, double> secondaries::NaivCompton::CalculateEnergy(
    double energy, double rho)
{
    auto energy_gamma = energy * (1 - rho);
    auto energy_electron = energy * rho;
    return make_tuple(energy_gamma, energy_electron);
}

vector<Loss::secondary_t> secondaries::NaivCompton::CalculateSecondaries(
    double primary_energy, Loss::secondary_t loss, const Component&,
    vector<double> rnd)
{
    auto rho = CalculateRho(primary_energy, get<Loss::ENERGY>(loss));
    auto secondary_energy = CalculateEnergy(get<Loss::ENERGY>(loss), rho);
    auto secondary_dir = CalculateDirections(
        get<Loss::DIRECTION>(loss), get<Loss::ENERGY>(loss), rho, rnd[1]);
    auto sec = std::vector<Loss::secondary_t>();
    sec.emplace_back(static_cast<int>(ParticleType::Gamma),
        get<Loss::POSITION>(loss), get<GAMMA>(secondary_dir),
        get<GAMMA>(secondary_energy), 0.);
    sec.emplace_back(static_cast<int>(ParticleType::Gamma),
        get<Loss::POSITION>(loss), get<ELECTRON>(secondary_dir),
        get<ELECTRON>(secondary_energy), 0.);
    return sec;
}
