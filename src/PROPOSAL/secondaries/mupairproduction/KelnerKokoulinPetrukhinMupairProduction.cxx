
#include "PROPOSAL/secondaries/mupairproduction/KelnerKokoulinPetrukhinMupairProduction.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using std::fmod;
using std::get;
using std::make_tuple;
using std::sqrt;

using namespace PROPOSAL;

secondaries::KelnerKokoulinPetrukhinMupairProduction::
    KelnerKokoulinPetrukhinMupairProduction(const ParticleDef& p, const Medium&)
    : p_def(p)
{
}

double secondaries::KelnerKokoulinPetrukhinMupairProduction::CalculateRho(
    double energy, double v, const Component& comp, double rnd)
{
    auto rho_max = 1 - 2 * MMU / (v * energy);
    if (rho_max < 0)
        return 0;
    integral.IntegrateWithRandomRatio(0, rho_max,
        [&](double rho) {
            return param.FunctionToIntegral(p_def, comp, energy, v, rho);
        },
        3, rnd);
    return integral.GetUpperLimit();
}

tuple<Vector3D, Vector3D>
secondaries::KelnerKokoulinPetrukhinMupairProduction::CalculateDirections(
    Vector3D primary_dir, double energy, double rho, double rnd)
{
    return make_tuple(primary_dir, primary_dir);
}

tuple<double, double>
secondaries::KelnerKokoulinPetrukhinMupairProduction::CalculateEnergy(
    double energy, double rho)
{
    auto energy_1 = 0.5 * energy * (1 + rho);
    auto energy_2 = 0.5 * energy * (1 - rho);
    return make_tuple(energy_1, energy_2);
}

vector<Loss::secondary_t>
secondaries::KelnerKokoulinPetrukhinMupairProduction::CalculateSecondaries(
    double initial_energy, Loss::secondary_t loss, const Component& comp,
    vector<double> rnd)
{
    auto v = get<Loss::ENERGY>(loss) / initial_energy;
    auto rho = CalculateRho(initial_energy, v, comp, rnd[0]);
    auto secondary_energy = CalculateEnergy(get<Loss::ENERGY>(loss), rho);
    auto secondary_dir = CalculateDirections(
        get<Loss::DIRECTION>(loss), get<Loss::ENERGY>(loss), rho, rnd[1]);
    auto sec = std::vector<Loss::secondary_t>();

    std::get<Loss::TYPE>(loss) = p_def.particle_type;
    std::get<Loss::ENERGY>(loss)
            = initial_energy - std::get<Loss::ENERGY>(loss);
    sec.emplace_back(loss);

    sec.emplace_back(static_cast<int>(ParticleType::MuMinus),
        get<Loss::POSITION>(loss), get<0>(secondary_dir),
        get<0>(secondary_energy), 0.);
    sec.emplace_back(static_cast<int>(ParticleType::MuPlus),
        get<Loss::POSITION>(loss), get<1>(secondary_dir),
        get<1>(secondary_energy), 0.);
    return sec;
}
