#include "PROPOSAL/secondaries/parametrization/annihilation/SingleDifferentialAnnihilation.h"
#include "PROPOSAL/Constants.h"

#include <cmath>
#include <stdexcept>
#include <tuple>

using std::fmod;
using std::make_tuple;
using std::sqrt;
using std::get;

using namespace PROPOSAL;

double secondaries::SingleDifferentialAnnihilation::CalculateRho(
    double energy, double rnd, const Component& comp)
{
    if (!dndx)
        throw std::logic_error("dndx Interpolant for SingleDifferentialAnnihilation not defined.");
    for (auto& it : *dndx) {
        if (comp.GetHash() == it.first) {
            auto rate = std::get<1>(it.second)->Calculate(energy);
            return std::get<1>(it.second)->GetUpperLimit(energy, rnd * rate);
        }
    }
    std::ostringstream s;
    s << "Component (" << comp.GetName()
      << ") can not be found in the precalculated annihilation tables.";
    throw std::out_of_range(s.str());
}

std::tuple<Cartesian3D, Cartesian3D>
secondaries::SingleDifferentialAnnihilation::CalculateDirections(
    const Vector3D& primary_dir, double energy, double rho, double rnd)
{
    auto com_energy = energy + ME; // center of mass energy
    auto kin_energy = energy - ME; // kinetic energy
    auto cosphi0 = (com_energy * (1. - rho) - ME)
        / ((1. - rho) * sqrt(com_energy * kin_energy));
    auto cosphi1
        = (com_energy * rho - ME) / (rho * sqrt(com_energy * kin_energy));
    auto rnd_theta = rnd * 2. * PI;
    auto dir_1 = Cartesian3D(primary_dir);
    auto dir_2 = Cartesian3D(primary_dir);
    dir_1.deflect(cosphi0, rnd_theta);
    dir_2.deflect(cosphi1, fmod(rnd_theta + PI, 2. * PI));
    return make_tuple(dir_1, dir_2);
}

std::tuple<double, double>
secondaries::SingleDifferentialAnnihilation::CalculateEnergy(
    double energy, double rho)
{
    assert(rho >= 0);
    assert(rho <= 1);
    auto energy_1 = (energy + ME) * (1 - rho);
    auto energy_2 = (energy + ME) * rho;
    return make_tuple(energy_1, energy_2);
}

std::vector<ParticleState>
secondaries::SingleDifferentialAnnihilation::CalculateSecondaries(
        StochasticLoss loss, const Component& comp, std::vector<double> &rnd)
{
    auto rho = CalculateRho(loss.parent_particle_energy, rnd[0], comp);
    auto secondary_energies = CalculateEnergy(loss.parent_particle_energy, rho);
    auto secondary_dir = CalculateDirections(
            loss.direction, loss.parent_particle_energy, rho, rnd[1]);
    auto sec = std::vector<ParticleState>();
    sec.emplace_back(ParticleType::Gamma, loss.position, get<0>(secondary_dir),
            get<0>(secondary_energies), loss.time, 0.);
    sec.emplace_back(ParticleType::Gamma, loss.position, get<1>(secondary_dir),
            get<1>(secondary_energies), loss.time, 0.);
    return sec;
}
