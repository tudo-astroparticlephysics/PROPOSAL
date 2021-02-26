
#include "PROPOSAL/secondaries/parametrization/mupairproduction/KelnerKokoulinPetrukhinMupairProduction.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using std::fmod;
using std::get;
using std::make_tuple;
using std::tuple;
using std::vector;
using std::sqrt;

using namespace PROPOSAL;


double secondaries::KelnerKokoulinPetrukhinMupairProduction::CalculateRho(
    double energy, double v, const Component& comp, double rnd1, double rnd2)
{
    auto rho_max = 1 - 2 * MMU / (v * energy);
    if (rho_max < 0)
        return 0;
    integral.IntegrateWithRandomRatio(0, rho_max,
        [&](double rho) {
            return param.FunctionToIntegral(p_def, comp, energy, v, rho);
        },
        3, rnd1);
    auto rho_tmp = integral.GetUpperLimit();
    if (rnd2 < 0.5) {
        return -rho_tmp;
    } else {
        return rho_tmp;
    }
}

tuple<Cartesian3D, Cartesian3D>
secondaries::KelnerKokoulinPetrukhinMupairProduction::CalculateDirections(
    const Vector3D& primary_dir, double, double, double)
{
    return make_tuple(Cartesian3D(primary_dir), Cartesian3D(primary_dir));
}

tuple<double, double>
secondaries::KelnerKokoulinPetrukhinMupairProduction::CalculateEnergy(
    double energy, double rho)
{
    auto energy_1 = 0.5 * energy * (1 + rho);
    auto energy_2 = 0.5 * energy * (1 - rho);
    return make_tuple(energy_1, energy_2);
}

vector<ParticleState>
secondaries::KelnerKokoulinPetrukhinMupairProduction::CalculateSecondaries(
    StochasticLoss loss, const Component& comp, vector<double> &rnd)
{
    auto v = loss.energy / loss.parent_particle_energy;
    auto rho = CalculateRho(loss.parent_particle_energy, v, comp, rnd[0], rnd[1]);
    auto secondary_energies = CalculateEnergy(loss.energy, rho);
    auto secondary_dir = CalculateDirections(
        loss.direction, loss.energy, rho, rnd[2]);

    auto sec = std::vector<ParticleState>();

    sec.emplace_back(static_cast<ParticleType>(p_def.particle_type),
                     loss.position, loss.direction,
                     loss.parent_particle_energy - loss.energy,
                     loss.time, loss.propagated_distance);

    sec.emplace_back(ParticleType::MuMinus, loss.position, get<0>(secondary_dir),
                     get<0>(secondary_energies), loss.time, 0.);

    sec.emplace_back(ParticleType::MuPlus, loss.position, get<1>(secondary_dir),
                     get<1>(secondary_energies), loss.time, 0.);

    return sec;
}
