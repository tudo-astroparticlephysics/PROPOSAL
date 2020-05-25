
#include <cmath>

#include "PROPOSAL/crossection/parametrization/MupairProduction.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"

#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/Constants.h"

#define MUPAIR_PARAM_INTEGRAL_IMPL(param)                                      \
    Mupair##param::Mupair##param()                                             \
        : MupairProductionRhoIntegral()                                        \
    {                                                                          \
    }

using namespace PROPOSAL;

MupairProduction::MupairProduction()
    : Parametrization(InteractionType::MuPair, "mupair")
    , drho_integral_(IROMB, IMAXS, IPREC)
{
}

Parametrization::KinematicLimits MupairProduction::GetKinematicLimits(
    const ParticleDef& p_def, const Component& comp, double energy)
{
    auto vmin = 2 * MMU / energy;
    auto vmax = 1 - p_def.mass / energy;

    if (vmax < vmin)
        vmax = vmin;

    return KinematicLimits(vmin, vmax);
}

double MupairProduction::Calculaterho(const ParticleDef& p_def, const Component& comp,
    double energy, double v, double rnd1, double rnd2)
{
    double rho = 0;
    double rho_min = 0;
    double rho_max = 1 - 2 * MMU / (v * energy);

    if (rho_max < 0)
        return 0;

    static_cast<void>(drho_integral_.IntegrateWithRandomRatio(rho_min, rho_max,
        std::bind(&MupairProduction::FunctionToIntegral, this, p_def, comp, energy, v,
            std::placeholders::_1),
        3, rnd1));

    rho = drho_integral_.GetUpperLimit();

    if (rnd2 < 0.5)
        rho = -rho;

    return rho;
}

MupairProductionRhoIntegral::MupairProductionRhoIntegral()
    : MupairProduction()
    , integral_(IROMB, IMAXS, IPREC)
{
}

double MupairProductionRhoIntegral::DifferentialCrossSection(const ParticleDef& p_def, const Component& comp,
    double energy, double v)
{
    double rMax, aux;

    aux = 1 - 2 * MMU / (v * energy);

    if (aux > 0) {
        rMax = aux;
    } else {
        return 0;
    }

    return NA / comp.GetAtomicNum() * p_def.charge
        * p_def.charge
        * (integral_.Integrate(0, rMax,
              std::bind(&MupairProductionRhoIntegral::FunctionToIntegral, this,
                 p_def, comp,  energy, v, std::placeholders::_1),
              2));
}

MUPAIR_PARAM_INTEGRAL_IMPL(KelnerKokoulinPetrukhin)

double MupairKelnerKokoulinPetrukhin::FunctionToIntegral(const ParticleDef& p_def, const Component& comp,
    double energy, double v, double r)
{
    // Parametrization of Kelner/Kokoulin/Petrukhin
    // Physics of Atomic Nuclei, Vol. 63, No. 9, 2000, pp. 1603-1611. Translated
    // from Yadernaya Fizika, Vol. 63, 2000, pp. 1690-1698 Original Russian Text
    // Copyright 2000 by Kel'ner, Kokoulin, Petukhin DOI: 10.1134/1.1312894

    double aux, aux1, aux2, r2, rMax, Z3, xi, beta, A_pow, r_mu;
    double phi, U, U_max, X, Y;
    double medium_charge = comp.GetNucCharge();
    double atomic_weight = comp.GetAtomInMolecule();
    // double medium_log_constant =
    // components_[component_index_]->GetLogConstant();
    double medium_log_constant = 183; // According to the paper, B is set to 183

    r2 = r * r;
    rMax = 1 - 2 * MMU / (v * energy);
    Z3 = std::pow(medium_charge, -1. / 3);
    aux = (p_def.mass * v) / (2 * MMU);
    xi = aux * aux * (1 - r2) / (1 - v);
    beta = (v * v) / (2 * (1 - v));
    A_pow = std::pow(atomic_weight, -0.27);
    r_mu = RE * ME / MMU; // classical muon radius

    // Phi Calculation (18)
    aux = (2 + r2) * (1 + beta) + xi * (3 + r2);
    aux *= std::log(1 + 1. / xi);

    aux1 = (1 + r2) * (1 + 1.5 * beta) - 1. / xi * (1 + 2 * beta) * (1 - r2);
    aux1 *= std::log(1 + xi);
    aux2 = -1 - 3 * r2 + beta * (1 - 2 * r2);

    phi = aux + aux1 + aux2;

    // X Calculation (22)
    Y = 12 * std::sqrt(MMU / energy); //(21)
    aux = 0.65 * A_pow * medium_log_constant * Z3 * MMU / ME;
    aux1 = 2 * SQRTE * std::pow(MMU, 2) * medium_log_constant * Z3 * (1 + xi)
        * (1 + Y);
    aux2 = ME * energy * v * (1 - r2);

    U = aux / (1 + aux1 / aux2);

    aux2 = ME * energy * v * (1 - rMax * rMax);
    U_max = aux / (1 + aux1 / aux2);

    X = 1 + U - U_max;

    // Combine results
    aux = ALPHA * r_mu * p_def.charge * medium_charge;
    aux *= 2 * aux * phi * (1 - v)
        / (1.5 * PI * v); // Factor 2: Similar to factor 2 from EPairProduction,
                          // probably from symmetry in Rho

    if (X > 0) {
        aux *= std::log(X);
    } else {
        aux = 0;
    }

    if (aux < 0) {
        return 0;
    }

    return aux;
}

#undef MUPAIR_PARAM_INTEGRAL_IMPL
