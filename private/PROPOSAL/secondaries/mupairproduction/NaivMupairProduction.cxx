
#include "PROPOSAL/secondaries/mupairproduction/NaivMupairProduction.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>
#include <stdexcept>

using std::fmod;
using std::make_tuple;
using std::sqrt;
using std::get;

using namespace PROPOSAL;

double secondaries::NaivMupairProduction::CalculateRho(
    double energy, double v, double rnd)
{
    auto rho_max = 1 - 2 * MMU / (v * energy);
    if (rho_max < 0)
        return 0;
    integral.IntegrateWithRandomRatio(0, rho_max,
        [&, energy, v](double rho) { return rho_integrand(energy, v, rho); }, 3,
        rnd);
    return integral.GetUpperLimit();
}

tuple<Vector3D, Vector3D>
secondaries::NaivMupairProduction::CalculateDirections(
    Vector3D primary_dir, double energy, double rho, double rnd)
{
    return make_tuple(primary_dir, primary_dir);
}

tuple<double, double> secondaries::NaivMupairProduction::CalculateEnergy(
    double energy, double rho)
{
    auto energy_1 = 0.5 * energy * (1 + rho);
    auto energy_2 = 0.5 * energy * (1 - rho);
    return make_tuple(energy_1, energy_2);
}

vector<Loss::secondary_t>
secondaries::NaivMupairProduction::CalculateSecondaries(Loss::secondary_t loss,
    array<double, secondaries::NaivMupairProduction::n_rnd> rnd)
{
    auto v = 0; // TODO: v initialization still missing;
    auto rho = CalculateRho(get<Loss::ENERGY>(loss), v, rnd[0]);
    auto secondary_energy = CalculateEnergy(get<Loss::ENERGY>(loss), rho);
    auto secondary_dir = CalculateDirections(
        get<Loss::DIRECTION>(loss), get<Loss::ENERGY>(loss), rho, rnd[1]);
    auto sec = std::vector<Loss::secondary_t>();
    sec.emplace_back(static_cast<int>(ParticleType::MuMinus),
        get<Loss::POSITION>(loss), get<0>(secondary_dir),
        get<0>(secondary_energy), 0.);
    sec.emplace_back(static_cast<int>(ParticleType::MuPlus),
        get<Loss::POSITION>(loss), get<1>(secondary_dir),
        get<1>(secondary_energy), 0.);
    return sec;
}
