#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/ionization/NaivIonization.h"

#include <cmath>
#include <stdexcept>

using std::fmod;
using std::make_tuple;
using std::sqrt;

using namespace PROPOSAL;

secondaries::NaivIonization::NaivIonization(const ParticleDef& p_def)
    : primary_particle_type(p_def.particle_type)
{
}

double secondaries::NaivIonization::CalculateRho(
    double primary_energy, double loss_energy)
{
    return primary_energy / loss_energy;
}

tuple<Vector3D, Vector3D> secondaries::NaivIonization::CalculateDirections(
    Vector3D, double, double, double)
{
    return make_tuple(dir, dir);
}

tuple<double, double> secondaries::NaivIonization::CalculateEnergy(
    double energy, double rho)
{
    return make_tuple(energy * (1 - rho), energy * rho);
}

virtual vector<Loss::secondary_t>
secondaries::NaivIonization::CalculateSecondaries(
    Loss::secondary_t, array<double, n_rnd>)
{
    auto rho = CalculateRho(get<Loss::ENERGY>(loss), rnd[0]);
    auto secondary_energy = CalculateEnergy(get<Loss::ENERGY>(loss), rho);
    auto secondary_dir = CalculateDirections(
        get<Loss::DIRECTION>(loss), get<Loss::ENERGY>(loss), rho, rnd[1]);
    auto sec = std::vector<Loss::secondary_t>();
    sec.emplace_back(primary_particle_type, get<Loss::POSITION>(loss),
        get<0>(secondary_dir), get<0>(secondary_energy), 0.);
    sec.emplace_back(static_cast<int>(ParticleType::Electron),
        get<Loss::POSITION>(loss), get<1>(secondary_dir),
        get<1>(secondary_energy), 0.);
    return sec;
}
