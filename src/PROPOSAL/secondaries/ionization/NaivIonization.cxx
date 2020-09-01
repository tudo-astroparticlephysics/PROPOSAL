#include "PROPOSAL/secondaries/ionization/NaivIonization.h"
#include "PROPOSAL/Constants.h"

#include <cmath>
#include <stdexcept>

using std::fmod;
using std::make_tuple;
using std::sqrt;

using namespace PROPOSAL;


double secondaries::NaivIonization::CalculateRho(
    double primary_energy, double loss_energy)
{
    return loss_energy / primary_energy;
}

tuple<Vector3D, Vector3D> secondaries::NaivIonization::CalculateDirections(
    Vector3D dir, double, double, double)
{
    return make_tuple(dir, dir);
}

tuple<double, double> secondaries::NaivIonization::CalculateEnergy(
    double energy, double rho)
{
    return make_tuple(energy * (1 - rho), energy * rho);
}

vector<Loss::secondary_t> secondaries::NaivIonization::CalculateSecondaries(
    double primary_energy, Loss::secondary_t loss, const Component&,
    vector<double> rnd)
{
    auto rho = CalculateRho(primary_energy, get<Loss::ENERGY>(loss));
    auto secondary_energy = CalculateEnergy(get<Loss::ENERGY>(loss), rho);
    auto secondary_dir = CalculateDirections(
        get<Loss::DIRECTION>(loss), get<Loss::ENERGY>(loss), rho, rnd[1]);
    auto sec = std::vector<Loss::secondary_t>();
    sec.emplace_back(primary_particle_type, get<Loss::POSITION>(loss),
        get<0>(secondary_dir), get<0>(secondary_energy), 0.);
    sec.emplace_back(static_cast<int>(ParticleType::EMinus),
        get<Loss::POSITION>(loss), get<1>(secondary_dir),
        get<1>(secondary_energy), 0.);
    return sec;
}
