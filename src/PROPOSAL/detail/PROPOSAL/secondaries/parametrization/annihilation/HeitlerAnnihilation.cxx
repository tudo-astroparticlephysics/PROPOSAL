#include "PROPOSAL/secondaries/parametrization/annihilation/HeitlerAnnihilation.h"
#include "PROPOSAL/math/MathMethods.h"

using namespace PROPOSAL;

double secondaries::HeitlerAnnihilation::f_term(double a1, double a2, double rho) {
    return a1 * std::log(rho) + a2 / rho - rho;
}

double secondaries::HeitlerAnnihilation::f(double rho, double a1, double a2, double rho_min, double rho_max, double rnd) {
    return f_term(a1, a2, rho) - f_term(a1, a2, rho_min) - rnd * (f_term(a1, a2, rho_max) - f_term(a1, a2, rho_min));
}

double secondaries::HeitlerAnnihilation::CalculateRho(double energy, double rnd, const Component&) {
    // solve the following equation for \rho:
    // \int_{\rho_{min}}^{\rho} \frac{d\sigma}{d\rho} = rnd * \int_{\rho_{min}}^{\rho_{max}} \frac{d\sigma}{d\rho}

    double gamma = energy / mass;
    double rho_min, rho_max;
    if (gamma <= 1) {
        rho_min = 0;
        rho_max = 0;
    } else {
        auto aux = std::sqrt((gamma - 1.) / (gamma + 1.));
        rho_min = 0.5 * (1. - aux);
        rho_max = 0.5 * (1. + aux);
    }

    double a2 = 1 / std::pow(gamma + 1, 2.);
    double a1 = 1 + 2 * gamma * a2;

    auto function_to_solve = [&](double rho) { return f(rho, a1, a2, rho_min, rho_max, rnd); };
    double precision = rho_min * energy * HALF_PRECISION;
    return Bisection(function_to_solve, rho_min, rho_max, precision, 100).first;
}

std::tuple<Cartesian3D, Cartesian3D>
secondaries::HeitlerAnnihilation::CalculateDirections(
        const Vector3D& primary_dir, double energy, double rho, double rnd)
{
    // use two-body interaction kinematics to infer angles of secondary particles
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
    return std::make_tuple(dir_1, dir_2);
}

std::tuple<double, double>
secondaries::HeitlerAnnihilation::CalculateEnergy(
        double energy, double rho)
{
    assert(rho >= 0);
    assert(rho <= 1);
    auto energy_1 = (energy + ME) * (1 - rho);
    auto energy_2 = (energy + ME) * rho;
    return std::make_tuple(energy_1, energy_2);
}

std::vector<ParticleState>
secondaries::HeitlerAnnihilation::CalculateSecondaries(
        StochasticLoss loss, const Component& comp, std::vector<double> &rnd)
{
    auto rho = CalculateRho(loss.parent_particle_energy, rnd[0], comp);
    auto secondary_energies = CalculateEnergy(loss.parent_particle_energy, rho);
    auto secondary_dir = CalculateDirections(
            loss.direction, loss.parent_particle_energy, rho, rnd[1]);
    auto sec = std::vector<ParticleState>();
    sec.emplace_back(ParticleType::Gamma, loss.position, std::get<0>(secondary_dir),
                     std::get<0>(secondary_energies), loss.time, 0.);
    sec.emplace_back(ParticleType::Gamma, loss.position, std::get<1>(secondary_dir),
                     std::get<1>(secondary_energies), loss.time, 0.);
    return sec;
}
