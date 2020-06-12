#include "PROPOSAL/secondaries/annihilation/NaivAnnihilation.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>
#include <stdexcept>

using std::fmod;
using std::make_tuple;
using std::sqrt;

using namespace PROPOSAL;

secondaries::NaivAnnihilation::NaivAnnihilation(
    const crosssection::Annihilation&, const ParticleDef&, const Medium&,
    const EnergyCutSettings&)
{
    throw std::logic_error("dNdx interpolants must be build.");
}

double secondaries::NaivAnnihilation::CalculateRho(double energy, double rnd)
{
    throw std::logic_error("Here must take a alogrithm take place which can "
                           "calculate out of given rnd the component "
                           "interaction take place and sample a random v.");
}

tuple<Vector3D, Vector3D> secondaries::NaivAnnihilation::CalculateDirections(
    Vector3D primary_dir, double energy, double rho, double rnd)
{
    auto com_energy = energy + ME; // center of mass energy
    auto kin_energy = energy - ME; // kinetic energy
    auto cosphi0 = (com_energy * (1. - rho) - ME)
        / ((1. - rho) * sqrt(com_energy * kin_energy));
    auto cosphi1
        = (com_energy * rho - ME) / (rho * sqrt(com_energy * kin_energy));
    auto rnd_theta = rnd * 2. * PI;
    auto dir_1 = deflect(primary_dir, cosphi0, rnd_theta);
    auto dir_2 = deflect(primary_dir, cosphi1, fmod(rnd_theta + PI, 2. * PI));
    return make_tuple(dir_1, dir_2);
}

tuple<double, double> secondaries::NaivAnnihilation::CalculateEnergy(
    double energy, double rho)
{
    auto energy_1 = (energy + ME) * (1 - rho);
    auto energy_2 = (energy + ME) * rho;
    return make_tuple(energy_1, energy_2);
}

vector<Loss::secondary_t> secondaries::NaivAnnihilation::CalculateSecondaries(
    Loss::secondary_t loss, array<double, secondaries::NaivAnnihilation::n_rnd> rnd)
{
    auto rho = CalculateRho(get<Loss::ENERGY>(loss), rnd[0]);
    auto secondary_energy = CalculateEnergy(get<Loss::ENERGY>(loss), rho);
    auto secondary_dir = CalculateDirections(
        get<Loss::DIRECTION>(loss), get<Loss::ENERGY>(loss), rho, rnd[1]);
    auto sec = std::vector<Loss::secondary_t>();
    sec.emplace_back(static_cast<int>(ParticleType::Gamma),
        get<Loss::POSITION>(loss), get<0>(secondary_dir),
        get<0>(secondary_energy), 0.);
    sec.emplace_back(static_cast<int>(ParticleType::Gamma),
        get<Loss::POSITION>(loss), get<1>(secondary_dir),
        get<1>(secondary_energy), 0.);
    return sec;
}
