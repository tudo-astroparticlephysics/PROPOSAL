
#include "PROPOSAL/secondaries/epairproduction/KelnerKokoulinPetrukhinEpairProduction.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>
#include <iostream>

using std::fmod;
using std::get;
using std::make_tuple;
using std::sqrt;

using namespace PROPOSAL;


double secondaries::KelnerKokoulinPetrukhinEpairProduction::CalculateRho(
        double energy, double v, const Component& comp, double rnd)
{

    auto aux = 1 - (4 * ME) / (energy * v);
    auto aux2 = 1 - (6 * p_def.mass * p_def.mass) / (energy * energy * (1 - v));

    double rho_max;
    if (aux > 0 && aux2 > 0) {
        rho_max = std::sqrt(aux) * aux2;
    } else {
        rho_max = 0;
    }

    integral.IntegrateWithRandomRatio(0, rho_max,
                                      [&](double rho) {
                                          return param.FunctionToIntegral(p_def, comp, energy, v, rho);
                                      },
                                      3, rnd);
    return integral.GetUpperLimit();
}

tuple<Vector3D, Vector3D>
secondaries::KelnerKokoulinPetrukhinEpairProduction::CalculateDirections(
        Vector3D primary_dir, double, double, double)
{
    return make_tuple(primary_dir, primary_dir);
}

tuple<double, double>
secondaries::KelnerKokoulinPetrukhinEpairProduction::CalculateEnergy(
        double energy, double rho)
{
    auto energy_1 = 0.5 * energy * (1 + rho);
    auto energy_2 = 0.5 * energy * (1 - rho);
    return make_tuple(energy_1, energy_2);
}

vector<Loss::secondary_t>
secondaries::KelnerKokoulinPetrukhinEpairProduction::CalculateSecondaries(
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

    sec.emplace_back(static_cast<int>(ParticleType::EMinus),
                     get<Loss::POSITION>(loss), get<0>(secondary_dir),
                     get<0>(secondary_energy), 0.);
    sec.emplace_back(static_cast<int>(ParticleType::EPlus),
                     get<Loss::POSITION>(loss), get<1>(secondary_dir),
                     get<1>(secondary_energy), 0.);
    return sec;
}
