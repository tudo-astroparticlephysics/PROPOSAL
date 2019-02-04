
#include <cmath>

#include "PROPOSAL/crossection/parametrization/MupairProduction.h"

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"

#define MUPAIR_PARAM_INTEGRAL_IMPL(param)                                                                               \
    Mupair##param::Mupair##param(const ParticleDef& particle_def,                                                        \
                               const Medium& medium,                                                                   \
                               const EnergyCutSettings& cuts,                                                          \
                               double multiplier,                                                                      \
                               bool mupair_enable)                                                                               \
        : MupairProductionRhoIntegral(particle_def, medium, cuts, multiplier, mupair_enable)                                      \
    {                                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    Mupair##param::Mupair##param(const Mupair##param& photo)                                                              \
        : MupairProductionRhoIntegral(photo)                                                                            \
    {                                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    Mupair##param::~Mupair##param() {}                                                                                   \
                                                                                                                       \
    const std::string Mupair##param::name_ = "Mupair" #param;

using namespace PROPOSAL;

/******************************************************************************
 *                              MupairProduction                               *
 ******************************************************************************/

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

MupairProduction::MupairProduction(const ParticleDef& particle_def,
                                 const Medium& medium,
                                 const EnergyCutSettings& cuts,
                                 double multiplier,
                                 bool mupair_enable)
    : Parametrization(particle_def, medium, cuts, multiplier)
    , mupair_enable_(mupair_enable)
{
}

MupairProduction::MupairProduction(const MupairProduction& mupair)
    : Parametrization(mupair)
    , mupair_enable_(mupair.mupair_enable_)
{
}

MupairProduction::~MupairProduction() {}

bool MupairProduction::compare(const Parametrization& parametrization) const
{
    const MupairProduction* pairproduction = static_cast<const MupairProduction*>(&parametrization);

    if (mupair_enable_ != pairproduction->mupair_enable_)
        return false;
    else
        return Parametrization::compare(parametrization);
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
Parametrization::IntegralLimits MupairProduction::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    limits.vMin = 2 * MMU / energy;
    limits.vMax = 1 - particle_def_.mass / energy;

    if (limits.vMax < limits.vMin)
    {
        limits.vMax = limits.vMin;
    }

    limits.vUp = std::min(limits.vMax, cut_settings_.GetCut(energy));

    if (limits.vUp < limits.vMin)
    {
        limits.vUp = limits.vMin;
    }

    return limits;
}


// ------------------------------------------------------------------------- //
// Parametrization of Kelner/Kokoulin/Petrukhin
// Proc. 12th ICCR (1971), 2436
// ------------------------------------------------------------------------- //

/******************************************************************************
 *                          MupairProduction Integral                           *
 ******************************************************************************/

// ------------------------------------------------------------------------- //
MupairProductionRhoIntegral::MupairProductionRhoIntegral(const ParticleDef& particle_def,
                                                       const Medium& medium,
                                                       const EnergyCutSettings& cuts,
                                                       double multiplier,
                                                       bool mupair_enable)
    : MupairProduction(particle_def, medium, cuts, multiplier, mupair_enable)
    , integral_(IROMB, IMAXS, IPREC)
{
}

// ------------------------------------------------------------------------- //
MupairProductionRhoIntegral::MupairProductionRhoIntegral(const MupairProductionRhoIntegral& mupair)
    : MupairProduction(mupair)
    , integral_(mupair.integral_)
{
}

// ------------------------------------------------------------------------- //
MupairProductionRhoIntegral::~MupairProductionRhoIntegral() {}

// ------------------------------------------------------------------------- //
bool MupairProductionRhoIntegral::compare(const Parametrization& parametrization) const
{
    const MupairProductionRhoIntegral* mupair = static_cast<const MupairProductionRhoIntegral*>(&parametrization);

    if (integral_ != mupair->integral_)
        return false;
    else
        return MupairProduction::compare(parametrization);
}

// ------------------------------------------------------------------------- //
double MupairProductionRhoIntegral::DifferentialCrossSection(double energy, double v)
{
    double rMax, aux;

    aux  = 1 - 2 * MMU / (v * energy);

    if (aux > 0)
    {
        rMax = aux;
    } else
    {
        rMax = 0;
    }

    aux = std::max(1 - rMax, COMPUTER_PRECISION);

    return multiplier_ * medium_->GetMolDensity() * components_[component_index_]->GetAtomInMolecule() *
           particle_def_.charge * particle_def_.charge *
           (integral_.Integrate(
                1 - rMax, aux, std::bind(&MupairProductionRhoIntegral::FunctionToIntegral, this, energy, v, std::placeholders::_1), 2) +
            integral_.Integrate(
                aux, 1, std::bind(&MupairProductionRhoIntegral::FunctionToIntegral, this, energy, v, std::placeholders::_1), 4));
}

// ------------------------------------------------------------------------- //
void MupairProductionRhoIntegral::print(std::ostream& os) const
{
    os << "mupair production enabled: " << mupair_enable_ << '\n';
}

// ------------------------------------------------------------------------- //
size_t MupairProductionRhoIntegral::GetHash() const
{
    size_t seed = Parametrization::GetHash();
    hash_combine(seed, mupair_enable_);

    return seed;
}

/******************************************************************************
 *                          Specifc Parametrizations                           *
 ******************************************************************************/

// ------------------------------------------------------------------------- //
// Define the specific parametrizations
// ------------------------------------------------------------------------- //

MUPAIR_PARAM_INTEGRAL_IMPL(KelnerKokoulinPetrukhin)

// ------------------------------------------------------------------------- //
double MupairKelnerKokoulinPetrukhin::FunctionToIntegral(double energy, double v, double r)
{
    // Parametrization of Kelner/Kokoulin/Petrukhin
    // Physics of Atomic Nuclei, Vol. 63, No. 9, 2000, pp. 1603-1611. Translated from Yadernaya Fizika, Vol. 63, 2000, pp. 1690-1698
    // Original Russian Text Copyright 2000 by Kel'ner, Kokoulin, Petukhin
    // DOI: 10.1134/1.1312894


    double g1, g2;
    double aux, aux1, aux2, r2;
    double diagram_e, diagram_mu, atomic_electron_contribution;
    double medium_charge       = components_[component_index_]->GetNucCharge();
    double medium_log_constant = components_[component_index_]->GetLogConstant();

    r           = 1 - r; // only for integral optimization - do not forget to swap integration limits!
    r2          = r * r;
    double Z3   = std::pow(medium_charge, -1. / 3);
    aux         = (particle_def_.mass * v) / (2 * ME);
    double xi   = aux * aux * (1 - r2) / (1 - v);
    double beta = (v * v) / (2 * (1 - v));

    // these are the Y_e and Y_mu expressions in the original paper
    diagram_e  = (5 - r2 + 4 * beta * (1 + r2)) / (2 * (1 + 3 * beta) * std::log(3 + 1 / xi) - r2 - 2 * beta * (2 - r2));
    diagram_mu = (4 + r2 + 3 * beta * (1 + r2)) / ((1 + r2) * (1.5 + 2 * beta) * std::log(3 + xi) + 1 - 1.5 * r2);

    // these arae the L_e and L_mu expressions
    aux        = (1.5 * ME) / (particle_def_.mass * Z3);
    aux1       = (1 + xi) * (1 + diagram_e);
    aux2       = (2 * ME * SQRTE * medium_log_constant * Z3) / (energy * v * (1 - r2));
    diagram_e  = std::log((medium_log_constant * Z3 * std::sqrt(aux1)) / (1 + aux2 * aux1)) - 0.5 * std::log(1 + aux * aux * aux1);
    diagram_mu = std::log((medium_log_constant / aux * Z3) / (1 + aux2 * (1 + xi) * (1 + diagram_mu)));

    // these are the Phi_e and Phi_mu expressions
    // if the logarithms above are below zero, the contribution of the diagram is set to zero
    // if lpm supression is taken into account, Phi_e is changed
    if (diagram_e > 0)
    {
        if (1 / xi < HALF_PRECISION)
        {
            // TODO: where does this short expression come from?
            diagram_e = (1.5 - r2 / 2 + beta * (1 + r2)) / xi * diagram_e;
        } else
        {
            diagram_e =
                (((2 + r2) * (1 + beta) + xi * (3 + r2)) * std::log(1 + 1 / xi) + (1 - r2 - beta) / (1 + xi) - (3 + r2)) *
                diagram_e;
        }
    } else
    {
        diagram_e = 0;
    }

    if (diagram_mu > 0)
    {
        diagram_mu = (((1 + r2) * (1 + 1.5 * beta) - (1 + 2 * beta) * (1 - r2) / xi) * std::log(1 + xi) +
                      xi * (1 - r2 - beta) / (1 + xi) + (1 + 2 * beta) * (1 - r2)) *
                     diagram_mu;
    } else
    {
        diagram_mu = 0;
    }

    // Calculating the contribution of atomic electrons as scattering partner
    // Phys. Atom. Nucl. 61 (1998), 448
    if (medium_charge == 1)
    {
        g1 = 4.4e-5;
        g2 = 4.8e-5;
    } else
    {
        g1 = 1.95e-5;
        g2 = 5.3e-5;
    }

    aux  = energy / particle_def_.mass;
    aux1 = 0.073 * std::log(aux / (1 + g1 * aux / (Z3 * Z3))) - 0.26;
    aux2 = 0.058 * std::log(aux / (1 + g2 * aux / Z3)) - 0.14;

    if (aux1 > 0 && aux2 > 0)
    {
        atomic_electron_contribution = aux1 / aux2;
    } else
    {
        atomic_electron_contribution = 0;
    }

    // combining the results
    aux = ALPHA * RE * particle_def_.charge;
    ;
    aux *= aux / (1.5 * PI) * 2 * medium_charge * (medium_charge + atomic_electron_contribution);
    aux1 = ME / particle_def_.mass * particle_def_.charge;
    ;
    aux *= (1 - v) / v * (diagram_e + aux1 * aux1 * diagram_mu);

    /*
    if (lpm_)
    {
        aux *= lpm(energy, v, r2, beta, xi);
    }
    */
    
    if (aux < 0)
    {
        aux = 0;
    }

    return aux;
}


#undef MUPAIR_PARAM_INTEGRAL_IMPL
