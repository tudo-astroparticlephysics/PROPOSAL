#include "PROPOSAL/secondaries/parametrization/photomupairproduction/PhotoMuPairProductionBurkhardtKelnerKokoulin.h"

using namespace PROPOSAL;

double secondaries::PhotoMuPairProductionBurkhardtKelnerKokoulin::Calculatex(
        double energy, double rnd, const Component &comp) {
    if (!dndx)
        throw std::logic_error("dndx Interpolant for PhotoMuPairProductionBurkhardtKelnerKokoulinnot defined.");
    for (auto& it : *dndx) {
        if (comp.GetHash() == it.first) {
            auto& calc = *std::get<1>(it.second);
            auto lim = calc.GetIntegrationLimits(energy);
            auto rate = rnd * calc.Calculate(energy, lim.max);
            auto rho = calc.GetUpperLimit(energy, rate);
            return rho;
        }
    }
    std::ostringstream s;
    s << "Component (" << comp.GetName()
    << ") can not be found in the precalculated photomupairproduction "
       "tables.";
    throw std::out_of_range(s.str());
}

std::tuple<Cartesian3D, Cartesian3D>
secondaries::PhotoMuPairProductionBurkhardtKelnerKokoulin::CalculateDirections(
        const Vector3D &dir, double energy, const Component &comp, double rnd1) {
    auto k = energy / MMU;

    auto cosphi0 = std::cos(1. / k);
    auto cosphi1 = std::cos(1. / k);

    auto theta0 = rnd1 * 2. * PI;
    auto theta1 = std::fmod(theta0 + PI, 2. * PI);
    if (cosphi0 == -1.)
        cosphi0 *= (-1);
    if (cosphi1 == -1.)
        cosphi1 *= (-1);
    auto dir_0 = Cartesian3D(dir);
    dir_0.deflect(cosphi0, theta0);
    auto dir_1 = Cartesian3D(dir);
    dir_1.deflect(cosphi1, theta1);
    return std::make_tuple(dir_0, dir_1);
}

std::vector<ParticleState>
secondaries::PhotoMuPairProductionBurkhardtKelnerKokoulin::CalculateSecondaries(
        StochasticLoss loss, const Component &comp, std::vector<double> &rnd) {
    auto x = Calculatex(loss.parent_particle_energy, rnd[0], comp);
    auto secondary_energies = CalculateEnergy(loss.parent_particle_energy, x);
    auto secondary_dir = CalculateDirections(
            loss.direction, loss.parent_particle_energy,
            comp, rnd[1]);

    auto sec = std::vector<ParticleState>();
    sec.emplace_back(ParticleType::MuMinus, loss.position,
                     std::get<0>(secondary_dir),
                     std::get<0>(secondary_energies), loss.time, 0.);
    sec.emplace_back(ParticleType::MuPlus, loss.position,
                     std::get<1>(secondary_dir),
                     std::get<1>(secondary_energies), loss.time, 0.);
    return sec;
}