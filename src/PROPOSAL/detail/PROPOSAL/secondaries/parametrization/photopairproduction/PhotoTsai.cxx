
#include "PROPOSAL/secondaries/parametrization/photopairproduction/PhotoTsai.h"
#include "PROPOSAL/Constants.h"

#include <cmath>
#include <tuple>

using std::make_tuple;

using namespace PROPOSAL;

using std::get;

double secondaries::PhotoTsai::FunctionToIntegral(
    double energy, double x, double theta, const Component& comp)
{

    // Pair production and bremsstrahlung of chraged leptons, Yung-Su Tsai,
    // Review of Modern Physics, Vol. 46, No. 4, October 1974
    // see formula (3.5)

    double aux;
    double E = energy * x; // electron energy
    double l = E * E * theta * theta / (ME * ME);
    double Z = comp.GetNucCharge();
    double Z3 = std::pow(comp.GetNucCharge(), -1. / 3);
    double G2 = Z * Z + Z;
    double tminprimesqrt = (ME * ME * (1. + l)) / (2. * energy * x * (1. - x));

    double z = std::pow(Z / 137., 2.);
    double f = 1.202 * z - 1.0369 * std::pow(z, 2.)
        + 1.008 * std::pow(z, 3.) / (1. + z); // (3.3)

    double delta;
    double B;
    double X;
    double Xel, Xinel;

    if (Z < 2.5) {
        double eta;

        if (Z < 1.5) {
            eta = 1.; // for hydrogen
        } else {
            eta = 1.6875; // for helium
        }

        delta = ME * ME / (2. * energy * x * (1. - x)); // (3.20)
        B = 2. * ALPHA * ME * eta / tminprimesqrt;      // (3.21)
        Xel = 2. * std::log(ME / delta) - std::log(1. + B * B) + 1. / 6.
            - (4. / 3.) / (1. + B * B)
            + (1. / 6.) / std::pow(1. + B * B, 2.); // (3.18) with erratum
        Xel *= Z * Z;
        Xinel = 2. * std::log(ME / delta) - std::log(1. + B * B) + 11. / 6.
            - 4. / (B * B) * std::log(1. + B * B) + (4. / 3.) / (1. + B * B)
            - (1. / 6.) / std::pow(1. + B * B, 2.); // (3.19) with erratum
        Xinel *= Z;
    } else {
        // a and aprime according to Table (B.4.)
        double a, aprime;
        if (Z <= 1.5) {
            a = 122.8;
            aprime = 282.4;
        } else if (Z <= 2.5) {
            a = 90.8;
            aprime = 265.8;
        } else if (Z <= 3.5) {
            a = 100.;
            aprime = 418.6;
        } else if (Z <= 4.5) {
            a = 106;
            aprime = 571.4;
        } else {
            a = 111.7;
            aprime = 724.2;
        }

        a *= Z3 / ME;
        aprime *= Z3 * Z3 / ME;

        Xel = std::log(a * a * ME * ME * std::pow(1. + l, 2.)
                  / (a * a * tminprimesqrt * tminprimesqrt + 1.))
            - 1.;
        Xel *= Z * Z; // (3.44)
        Xinel = std::log(aprime * aprime * ME * ME * std::pow(1. + l, 2.)
                    / (aprime * aprime * tminprimesqrt * tminprimesqrt + 1.))
            - 1.;
        Xinel *= Z; // (3.45)
    }

    X = Xel + Xinel;

    // (3.5) with erratum
    aux = G2
        * (2. * x * (1. - x) / std::pow(1. + l, 2.)
            - 12. * l * x * (1. - x) / std::pow(1. + l, 4.));
    aux += (X - 2. * Z * Z * f)
        * ((2. * x * x - 2. * x + 1.) / std::pow(1. + l, 2.)
            + (4. * l * x * (1. - x)) / std::pow(1. + l, 4));

    // aux *= 2. * std::pow(ALPHA, 3.) / (PI * energy) * (E * E /
    // std::pow(ME, 4.)); only overall factor relevant

    aux *= sin(theta); // conversion from differential in cos(theta) to
                       // differential in theta

    return aux;
}

double secondaries::PhotoTsai::CalculateRho(
    double energy, double rnd, const Component& comp)
{
    if (!dndx)
        throw std::logic_error("dndx Interpolant for PhotoTsai not defined.");
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
      << ") can not be found in the precalculated tsai photopairproduction "
         "tables.";
    throw std::out_of_range(s.str());
}

std::tuple<Cartesian3D, Cartesian3D> secondaries::PhotoTsai::CalculateDirections(
    const Vector3D& dir, double energy, double rho, const Component& comp,
    std::vector<double> rnd)
{
    auto subst = std::max(1., std::log10(energy));
    auto integrand_substitution = [&, energy, rho, comp](double t) {
        return subst * std::pow(t, subst - 1.)
            * FunctionToIntegral(energy, rho, std::pow(t, subst), comp);
    };
    auto t_max = std::pow(PI, 1. / subst);
    integral.IntegrateWithRandomRatio(
        0., t_max, integrand_substitution, 3, rnd[2]);
    auto cosphi0 = std::cos(std::pow(integral.GetUpperLimit(), subst));
    integral.IntegrateWithRandomRatio(
        0., t_max, integrand_substitution, 3, rnd[3]);
    auto cosphi1 = std::cos(std::pow(integral.GetUpperLimit(), subst));
    auto theta0 = rnd[4] * 2. * PI;
    auto theta1 = std::fmod(theta0 + PI, 2. * PI);
    // TODO: Sometimes the intergration fails and -1 instead of 1 is returned...
    // :(
    if (cosphi0 == -1.)
        cosphi0 *= (-1);
    if (cosphi1 == -1.)
        cosphi1 *= (-1);
    auto dir_0 = Cartesian3D(dir);
    dir_0.deflect(cosphi0, theta0);
    auto dir_1 = Cartesian3D(dir);
    dir_1.deflect(cosphi1, theta1);
    return make_tuple(dir_0, dir_1);
}

std::tuple<double, double> secondaries::PhotoTsai::CalculateEnergy(
    double energy, double rho, double rnd)
{
    if (rnd > 0.5)
        return make_tuple(energy * (1 - rho), energy * rho);
    return make_tuple(energy * rho, energy * (1 - rho));
}

std::vector<ParticleState> secondaries::PhotoTsai::CalculateSecondaries(
    StochasticLoss loss, const Component& comp, std::vector<double>& rnd)
{
    auto rho = CalculateRho(loss.parent_particle_energy, rnd[0], comp);
    auto secondary_energies
        = CalculateEnergy(loss.parent_particle_energy, rho, rnd[1]);
    auto secondary_dir = CalculateDirections(
        loss.direction, loss.parent_particle_energy, rho, comp, rnd);

    auto sec = std::vector<ParticleState>();
    sec.emplace_back(ParticleType::EMinus, loss.position, get<0>(secondary_dir),
        get<0>(secondary_energies), loss.time, 0.);
    sec.emplace_back(ParticleType::EPlus, loss.position, get<1>(secondary_dir),
        get<1>(secondary_energies), loss.time, 0.);
    return sec;
}
