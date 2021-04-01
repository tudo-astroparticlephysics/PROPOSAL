
#include "PROPOSAL/secondaries/parametrization/epairproduction/KelnerKokoulinPetrukhinEpairProduction.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>
#include <sstream>

using std::fmod;
using std::get;
using std::make_tuple;
using std::sqrt;

using namespace PROPOSAL;

double secondaries::KelnerKokoulinPetrukhinEpairProduction::CalculateRho(
    double energy, double v, const Component& comp, double rnd1, double rnd2)
{

    auto aux = 1 - (4 * ME) / (energy * v);
    auto aux2 = 1 - (6 * p_def.mass * p_def.mass) / (energy * energy * (1 - v));

    double rho_max = std::sqrt(aux) * aux2;

    if (rho_max < 0) {
        std::stringstream ss;
        ss << "Rho should never be smaller than zero. Something with aux: "
           << aux << ", aux2: " << aux2 << " got wrong.";
        throw std::logic_error(ss.str());
    }

    auto integrand = [&](double rho) {
        return param.FunctionToIntegral(p_def, comp, energy, v, rho);
    };

    auto rho = 0.;
    if (integral.IntegrateWithRandomRatio(0, rho_max, integrand, 3, rnd1) > 0)
        rho = integral.GetUpperLimit();

    if (rnd2 < 0.5)
        rho *= -1.;

    return rho;
}

std::tuple<Cartesian3D, Cartesian3D>
secondaries::KelnerKokoulinPetrukhinEpairProduction::CalculateDirections(
    const Vector3D& primary_dir, double, double, double)
{
    return make_tuple(Cartesian3D(primary_dir), Cartesian3D(primary_dir));
}

std::tuple<double, double>
secondaries::KelnerKokoulinPetrukhinEpairProduction::CalculateEnergy(
    double energy, double rho)
{
    auto energy_1 = 0.5 * energy * (1 + rho);
    auto energy_2 = 0.5 * energy * (1 - rho);
    return make_tuple(energy_1, energy_2);
}

std::vector<ParticleState>
secondaries::KelnerKokoulinPetrukhinEpairProduction::CalculateSecondaries(
    StochasticLoss loss, const Component& comp, std::vector<double>& rnd)
{
    auto v = loss.energy / loss.parent_particle_energy;
    auto rho
        = CalculateRho(loss.parent_particle_energy, v, comp, rnd[0], rnd[1]);
    auto secondary_energies = CalculateEnergy(loss.energy, rho);
    auto secondary_dir
        = CalculateDirections(loss.direction, loss.energy, rho, rnd[2]);

    auto sec = std::vector<ParticleState>();
    sec.emplace_back(static_cast<ParticleType>(p_def.particle_type),
        loss.position, loss.direction,
        loss.parent_particle_energy - loss.energy, loss.time,
        loss.propagated_distance);

    sec.emplace_back(ParticleType::EMinus, loss.position, get<0>(secondary_dir),
        get<0>(secondary_energies), loss.time, 0.);

    sec.emplace_back(ParticleType::EPlus, loss.position, get<1>(secondary_dir),
        get<1>(secondary_energies), loss.time, 0.);

    return sec;
}
