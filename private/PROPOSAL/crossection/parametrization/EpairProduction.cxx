
#include <boost/bind.hpp>
#include <boost/functional/hash.hpp>
#include <cmath>

#include "PROPOSAL/crossection/parametrization/EpairProduction.h"

#include "PROPOSAL/math//Interpolant.h"
#include "PROPOSAL/math//InterpolantBuilder.h"

#include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"

using namespace PROPOSAL;

/******************************************************************************
*                              EpairProduction                               *
******************************************************************************/

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

EpairProduction::EpairProduction(const ParticleDef& particle_def,
                                 const Medium& medium,
                                 const EnergyCutSettings& cuts,
                                 double multiplier,
                                 bool lpm)
    : Parametrization(particle_def, medium, cuts, multiplier)
    , v_(0)
    , init_lpm_effect_(true)
    , lpm_(lpm)
    , eLpm_(0)
{
}

EpairProduction::EpairProduction(const EpairProduction& epair)
    : Parametrization(epair)
    , v_(epair.v_)
    , init_lpm_effect_(true)
    , lpm_(true)
    , eLpm_(epair.eLpm_)
{
}

EpairProduction::~EpairProduction()
{
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
Parametrization::IntegralLimits EpairProduction::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    double aux;
    double particle_energy = energy;
    double particle_mass   = particle_def_.mass;

    limits.vMin = 4 * ME / particle_energy;
    limits.vMax =
        1 - (3. / 4) * SQRTE * (particle_mass / particle_energy) * pow(components_[component_index_]->GetNucCharge(), 1. / 3);
    aux         = particle_mass / particle_energy;
    aux         = 1 - 6 * aux * aux;
    limits.vMax = std::min(limits.vMax, aux);
    limits.vMax = std::min(limits.vMax, 1 - particle_mass / particle_energy);

    if (limits.vMax < limits.vMin)
    {
        limits.vMax = limits.vMin;
    }

    limits.vUp = std::min(limits.vMax, cut_settings_.GetCut(particle_energy));

    if (limits.vUp < limits.vMin)
    {
        limits.vUp = limits.vMin;
    }

    return limits;
}

// ------------------------------------------------------------------------- //
double EpairProduction::FunctionToIntegral(double energy, double r)
{

    double Fe, Fm, Le, Lm, Ye, Ym, s, b, g1, g2;
    double aux, aux1, aux2, r2, Z3, atomic_electron_contribution;
    double particle_mass       = particle_def_.mass;
    double particle_charge     = particle_def_.charge;
    double particle_energy     = energy;
    double medium_charge       = components_[component_index_]->GetNucCharge();
    double medium_log_constant = components_[component_index_]->GetLogConstant();

    r   = 1 - r; // only for integral optimization - do not forget to swap integration limits!
    r2  = r * r;
    Z3  = pow(medium_charge, -1. / 3);
    aux = (particle_mass * v_) / (2 * ME);
    aux *= aux;
    s    = aux * (1 - r2) / (1 - v_);
    b    = (v_ * v_) / (2 * (1 - v_));
    Ye   = (5 - r2 + 4 * b * (1 + r2)) / (2 * (1 + 3 * b) * log(3 + 1 / s) - r2 - 2 * b * (2 - r2));
    Ym   = (4 + r2 + 3 * b * (1 + r2)) / ((1 + r2) * (1.5 + 2 * b) * log(3 + s) + 1 - 1.5 * r2);
    aux  = (1.5 * ME) / (particle_mass * Z3);
    aux1 = (1 + s) * (1 + Ye);
    aux2 = (2 * ME * SQRTE * medium_log_constant * Z3) / (particle_energy * v_ * (1 - r2));
    Le   = log((medium_log_constant * Z3 * sqrt(aux1)) / (1 + aux2 * aux1)) - 0.5 * log(1 + aux * aux * aux1);
    Lm   = log((medium_log_constant / aux * Z3) / (1 + aux2 * (1 + s) * (1 + Ym)));

    if (Le > 0)
    {
        if (1 / s < HALF_PRECISION)
        {
            // TODO: where does this short expression come from?
            Fe = (1.5 - r2 / 2 + b * (1 + r2)) / s * Le;
        } else
        {
            Fe = (((2 + r2) * (1 + b) + s * (3 + r2)) * log(1 + 1 / s) + (1 - r2 - b) / (1 + s) - (3 + r2)) * Le;
        }
    } else
    {
        Fe = 0;
    }

    if (Lm > 0)
    {
        Fm = (((1 + r2) * (1 + 1.5 * b) - (1 + 2 * b) * (1 - r2) / s) * log(1 + s) + s * (1 - r2 - b) / (1 + s) +
              (1 + 2 * b) * (1 - r2)) *
             Lm;
    }

    else
    {
        Fm = 0;
    }

    // Calculating the contribution of atomic electrons as scattering partner
    if (medium_charge == 1)
    {
        g1 = 4.4e-5;
        g2 = 4.8e-5;
    } else
    {
        g1 = 1.95e-5;
        g2 = 5.3e-5;
    }

    aux  = particle_energy / particle_mass;
    aux1 = 0.073 * log(aux / (1 + g1 * aux / (Z3 * Z3))) - 0.26;
    aux2 = 0.058 * log(aux / (1 + g2 * aux / Z3)) - 0.14;

    if (aux1 > 0 && aux2 > 0)
    {
        atomic_electron_contribution = aux1 / aux2;
    } else
    {
        atomic_electron_contribution = 0;
    }

    aux = ALPHA * RE;
    aux *= aux / (1.5 * PI);
    aux1 = ME / particle_mass * particle_charge;
    aux1 *= aux1;
    aux *= 2 * medium_charge * (medium_charge + atomic_electron_contribution) * ((1 - v_) / v_) * (Fe + aux1 * Fm);

    if (lpm_)
    {
        aux *= lpm(energy, r2, b, s);
    }

    if (aux < 0)
    {
        aux = 0;
    }

    return aux;
}

// ------------------------------------------------------------------------- //
double EpairProduction::lpm(double energy, double r2, double b, double x)
{
    if (init_lpm_effect_)
    {
        init_lpm_effect_ = false;
        double sum       = 0;

        for (int i = 0; i < medium_->GetNumComponents(); i++)
        {
            Components::Component* component = medium_->GetComponents().at(i);

            sum += component->GetNucCharge() * component->GetNucCharge() *
                   log(3.25 * component->GetLogConstant() * pow(component->GetNucCharge(), -1. / 3));
        }
        double particle_mass   = particle_def_.mass;
        double particle_charge = particle_def_.charge;

        eLpm_ = particle_mass / (ME * RE);
        eLpm_ *= (eLpm_ * eLpm_) * ALPHA * particle_mass /
                 (2 * PI * medium_->GetMolDensity() * particle_charge * particle_charge * sum);
    }

    double A, B, C, D, E, s;
    double s2, s36, s6, d1, d2, atan_, log1, log2;

    s     = sqrt(eLpm_ / (energy * v_ * (1 - r2))) / 4;
    s6    = 6 * s;
    atan_ = s6 * (x + 1);

    if (atan_ > 1 / COMPUTER_PRECISION)
    {
        return 1;
    }

    s2    = s * s;
    s36   = 36 * s2;
    d1    = s6 / (s6 + 1);
    d2    = s36 / (s36 + 1);
    atan_ = atan(atan_) - PI / 2;
    log1  = log((s36 * (1 + x) * (1 + x) + 1) / (s36 * x * x));
    log2  = log((s6 * (1 + x) + 1) / (s6 * x));
    A     = 0.5 * d2 * (1 + 2 * d2 * x) * log1 - d2 + 6 * d2 * s * (1 + ((s36 - 1) / (s36 + 1)) * x) * atan_;
    B     = d1 * (1 + d1 * x) * log2 - d1;
    C     = -d2 * d2 * x * log1 + d2 - (d2 * d2 * (s36 - 1) / (6 * s)) * x * atan_;
    D     = d1 - d1 * d1 * x * log2;
    E     = -s6 * atan_;

    return ((1 + b) * (A + (1 + r2) * B) + b * (C + (1 + r2) * D) + (1 - r2) * E) /
           (((2 + r2) * (1 + b) + x * (3 + r2)) * log(1 + 1 / x) + (1 - r2 - b) / (1 + x) - (3 + r2));
}

// ------------------------------------------------------------------------- //
// Getter
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
const std::string EpairProduction::name_ = "EpairProduction";

// ------------------------------------------------------------------------- //
size_t EpairProduction::GetHash() const
{
    size_t seed = Parametrization::GetHash();
    boost::hash_combine(seed, lpm_);

    return seed;
}

/******************************************************************************
*                         EpairProductionRhoIntegral                         *
******************************************************************************/


// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

EpairProductionRhoIntegral::EpairProductionRhoIntegral(const ParticleDef& particle_def,
                                 const Medium& medium,
                                 const EnergyCutSettings& cuts,
                                 double multiplier,
                                 bool lpm)
    : EpairProduction(particle_def, medium, cuts, multiplier, lpm)
    , integral_(IROMB, IMAXS, IPREC)
{
}

EpairProductionRhoIntegral::EpairProductionRhoIntegral(const EpairProductionRhoIntegral& epair)
    : EpairProduction(epair)
    , integral_(epair.integral_)
{
}

EpairProductionRhoIntegral::~EpairProductionRhoIntegral()
{
}

// ------------------------------------------------------------------------- //
// Crosssection
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double EpairProductionRhoIntegral::DifferentialCrossSection(double energy, double v)
{
    double particle_energy = energy;

    double rMax, aux, aux2;
    double particle_charge = particle_def_.charge;
    double particle_mass   = particle_def_.mass;

    v_   = v;
    aux  = 1 - (4 * ME) / (particle_energy * v_);
    aux2 = 1 - (6 * particle_mass * particle_mass) / (particle_energy * particle_energy * (1 - v_));

    if (aux > 0 && aux2 > 0)
    {
        rMax = sqrt(aux) * aux2;
    } else
    {
        rMax = 0;
    }

    aux = std::max(1 - rMax, COMPUTER_PRECISION);

    return multiplier_ * medium_->GetMolDensity() * components_[component_index_]->GetAtomInMolecule() * particle_charge * particle_charge *
           (integral_.Integrate(1 - rMax, aux, boost::bind(&EpairProductionRhoIntegral::FunctionToIntegral, this, energy, _1), 2) +
            integral_.Integrate(aux, 1, boost::bind(&EpairProductionRhoIntegral::FunctionToIntegral, this, energy, _1), 4));
}

/******************************************************************************
*                         EpairProductionRhoInterpolant                       *
******************************************************************************/

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

EpairProductionRhoInterpolant::EpairProductionRhoInterpolant(const ParticleDef& particle_def,
                                 const Medium& medium,
                                 const EnergyCutSettings& cuts,
                                 double multiplier,
                                 bool lpm,
                                 InterpolationDef def)
    : EpairProductionRhoIntegral(particle_def, medium, cuts, multiplier, lpm)
    , interpolant_(medium_->GetNumComponents(), NULL)
{
    std::vector<Interpolant2DBuilder> builder2d(components_.size());
    Helper::InterpolantBuilderContainer builder_container2d(components_.size());

    for (unsigned int i = 0; i < components_.size(); ++i)
    {
        builder2d[i].SetMax1(NUM1)
            .SetX1Min(particle_def_.low)
            .SetX1Max(BIGENERGY)
            .SetMax2(NUM1)
            .SetX2Min(0.0)
            .SetX2Max(1.0)
            .SetRomberg1(def.order_of_interpolation)
            .SetRational1(false)
            .SetRelative1(false)
            .SetIsLog1(true)
            .SetRomberg2(def.order_of_interpolation)
            .SetRational2(false)
            .SetRelative2(false)
            .SetIsLog2(false)
            .SetRombergY(def.order_of_interpolation)
            .SetRationalY(false)
            .SetRelativeY(false)
            .SetLogSubst(false)
            .SetFunction2D(boost::bind(&EpairProductionRhoInterpolant::FunctionToBuildEpairInterpolant, this, _1, _2, i));

        builder_container2d[i].first = &builder2d[i];
        builder_container2d[i].second = &interpolant_[i];
    }

    Helper::InitializeInterpolation("Epair",
                                    builder_container2d,
                                    std::vector<Parametrization*>(1, this),
                                    def);
}

EpairProductionRhoInterpolant::EpairProductionRhoInterpolant(const EpairProductionRhoInterpolant& epair)
    : EpairProductionRhoIntegral(epair)
    , interpolant_()
{
    interpolant_.resize(epair.interpolant_.size());

    for(unsigned int i = 0; i < epair.interpolant_.size(); ++i)
    {
        interpolant_[i] = new Interpolant(*epair.interpolant_[i]);
    }
}

EpairProductionRhoInterpolant::~EpairProductionRhoInterpolant()
{
    for(unsigned int i = 0; i < interpolant_.size(); ++i)
    {
        delete interpolant_[i];
    }

    interpolant_.clear();
}

// ------------------------------------------------------------------------- //
// Crosssection
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double EpairProductionRhoInterpolant::DifferentialCrossSection(double energy, double v)
{
    Parametrization::IntegralLimits limits = GetIntegralLimits(energy);

    if (v >= limits.vUp)
    {
        return std::max(
            interpolant_.at(component_index_)->Interpolate(energy, log(v / limits.vUp) / log(limits.vMax / limits.vUp)), 0.0);
    }
    else
    {
        return EpairProductionRhoIntegral::DifferentialCrossSection(energy, v);
    }
}

// ------------------------------------------------------------------------- //
double EpairProductionRhoInterpolant::FunctionToBuildEpairInterpolant(double energy, double v, int component)
{
    component_index_ = component;
    Parametrization::IntegralLimits limits = GetIntegralLimits(energy);

    if (limits.vUp == limits.vMax)
    {
        return 0;
    }

    v = limits.vUp * exp(v * log(limits.vMax / limits.vUp));

    return EpairProductionRhoIntegral::DifferentialCrossSection(energy, v);
}
