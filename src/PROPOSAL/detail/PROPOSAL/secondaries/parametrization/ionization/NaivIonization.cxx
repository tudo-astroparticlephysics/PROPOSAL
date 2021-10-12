#include "PROPOSAL/secondaries/parametrization/ionization/NaivIonization.h"
#include "PROPOSAL/Constants.h"

#include <cmath>
#include <stdexcept>

using std::fmod;
using std::make_tuple;
using std::sqrt;
using std::get;

using namespace PROPOSAL;



std::tuple<Cartesian3D, Cartesian3D> secondaries::NaivIonization::CalculateDirections(
    const Vector3D& dir, double E_i, double v, double rnd)
{
    // Calculate directions using momentum conservation
    auto E_f = E_i * (1 - v);
    auto E_e = E_i * v + ME;

    auto p_i = std::sqrt((E_i + primary_particle_mass) * (E_i - primary_particle_mass));
    auto p_f = std::sqrt((E_f + primary_particle_mass) * (E_f - primary_particle_mass));
    auto p_e = std::sqrt((E_e + ME) * (E_e - ME));

    // deflection angle between ingoing and outgoing primary particle
    auto cosphi = ((E_i + ME) * E_f - E_i * ME - primary_particle_mass * primary_particle_mass) / (p_i * p_f);
    auto dir_outgoing = Cartesian3D(dir);
    dir_outgoing.deflect(cosphi, rnd * 2. * PI);

    // deflection angle between ingoing particle and outgoing delta ray
    auto cosphi_bar = ((E_i + ME) * E_e - E_i * ME - ME * ME) / (p_i * p_e);
    auto dir_delta = Cartesian3D(dir);
    dir_delta.deflect(cosphi_bar, std::fmod(rnd * 2. * PI + PI, 2. * PI));

    return make_tuple(dir_outgoing, dir_delta);
}

std::tuple<double, double> secondaries::NaivIonization::CalculateEnergy(
    double energy, double v)
{
    return make_tuple(energy * (1 - v), energy * v + ME);
}

std::vector<ParticleState> secondaries::NaivIonization::CalculateSecondaries(
    StochasticLoss loss, const Component&, std::vector<double> &rnd)
{
    auto v = loss.energy / loss.parent_particle_energy;
    auto secondary_energies = CalculateEnergy(loss.parent_particle_energy, v);
    auto secondary_dir = CalculateDirections(
            loss.direction, loss.parent_particle_energy, v, rnd[0]);

    auto sec = std::vector<ParticleState>();
    sec.emplace_back(static_cast<ParticleType>(primary_particle_type),
                     loss.position, get<0>(secondary_dir),
                     get<0>(secondary_energies), loss.time, 0.);
    sec.emplace_back(ParticleType::EMinus, loss.position, get<1>(secondary_dir),
                     get<1>(secondary_energies), loss.time, 0.);
    return sec;
}
