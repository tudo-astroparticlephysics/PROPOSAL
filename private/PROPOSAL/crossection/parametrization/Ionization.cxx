
#include <boost/bind.hpp>
#include <cmath>

#include "PROPOSAL/crossection/parametrization/Ionization.h"

#include "PROPOSAL/Output.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/Integral.h"

using namespace PROPOSAL;

/******************************************************************************
*                               Ionization                                *
******************************************************************************/

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

Ionization::Ionization(const ParticleDef& particle_def,
                       const Medium& medium,
                       const EnergyCutSettings& cuts,
                       double multiplier)
    : Parametrization(particle_def, medium, cuts, multiplier)
{
}

Ionization::Ionization(const Ionization& ioniz)
    : Parametrization(ioniz)
{
}

Ionization::~Ionization()
{
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double Ionization::DifferentialCrossSection(double energy, double v)
{
    double result, aux, aux2;

    IntegralLimits limits = GetIntegralLimits(energy);

    // TODO(mario): Better way? Sat 2017/09/02
    double square_momentum   = energy * energy - particle_def_.mass * particle_def_.mass;
    double particle_momentum = sqrt(std::max(square_momentum, 0.0));
    double beta              = particle_momentum / energy;
    double gamma             = energy / particle_def_.mass;

    aux  = beta * beta;
    aux2 = v / (1 + 1 / gamma);
    aux2 *= 0.5 * aux2;
    result = 1 - aux * (v / limits.vMax) + aux2;
    result *= IONK * particle_def_.charge * particle_def_.charge * medium_->GetZA() / (2 * aux * energy * v * v);

    return multiplier_ * medium_->GetMassDensity() * result;
}

// ------------------------------------------------------------------------- //
double Ionization::FunctionToDEdxIntegral(double energy, double variable)
{
    return variable * DifferentialCrossSection(energy, variable) * InelCorrection(energy, variable);
}

//----------------------------------------------------------------------------//
double Ionization::FunctionToDNdxIntegral(double energy, double variable)
{
    return DifferentialCrossSection(energy, variable) * (1 + InelCorrection(energy, variable));
}

// ------------------------------------------------------------------------- //
Parametrization::IntegralLimits Ionization::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    double aux;
    double gamma = energy / particle_def_.mass;

    limits.vMin = (1.e-6 * medium_->GetI()) / energy;
    aux         = ME / particle_def_.mass;
    limits.vMax = 2 * ME * (gamma * gamma - 1) / ((1 + 2 * gamma * aux + aux * aux) * energy);
    limits.vMax = std::min(limits.vMax, 1. - particle_def_.mass / energy);

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
double Ionization::InelCorrection(double energy, double v)
{
    double result, a, b, c;
    IntegralLimits limits = GetIntegralLimits(energy);

    double gamma = energy / particle_def_.mass;

    a      = log(1 + 2 * v * energy / ME);
    b      = log((1 - v / limits.vMax) / (1 - v));
    c      = log((2 * gamma * (1 - v) * ME) / (particle_def_.mass * v));
    result = a * (2 * b + c) - b * b;

    return (ALPHA / (2 * PI)) * result;
}

// ------------------------------------------------------------------------- //
// Getter
// ------------------------------------------------------------------------- //

const std::string Ionization::name_ = "Ionization";
