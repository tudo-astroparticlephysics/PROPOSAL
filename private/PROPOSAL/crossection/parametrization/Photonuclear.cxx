
#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <boost/functional/hash.hpp>
#include <cmath>

#include "PROPOSAL/crossection/parametrization/Photonuclear.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/Components.h"


using namespace PROPOSAL;

/******************************************************************************
*                                 RealPhoton                                  *
******************************************************************************/

bool RealPhoton::operator==(const RealPhoton& photon) const
{
    if (typeid(*this) != typeid(photon))
        return false;
    else
        return compare(photon);
}

bool RealPhoton::operator!=(const RealPhoton& photon) const
{
    return !(*this == photon);
}

bool RealPhoton::compare(const RealPhoton& photon) const
{
    (void) photon;
    return true;
}


/******************************************************************************
*                                   HardBB                                    *
******************************************************************************/

std::vector<double> HardBB::x = boost::assign::list_of(3)(4)(5)(6)(7)(8)(9);

const std::string HardBB::name_ = "HardBB";
const std::string SoftBB::name_ = "SoftBB";

HardBB::HardBB(const ParticleDef& particle_def)
    : interpolant_()
{
    const HardBBTables::VecType& y = particle_def.hardbb_table;

    if (!y.empty())
    {
        for(unsigned int i=0; i < y.size(); i++)
        {
            interpolant_.push_back(new Interpolant(x, y.at(i), 4, false, false));
        }
    }
    else
    {
        log_fatal("No HardBB tables provided for the given particle %s", particle_def.name.c_str());
    }
}

HardBB::HardBB(const HardBB& hardBB)
    : interpolant_()
{
    interpolant_.resize(hardBB.interpolant_.size());

    for(unsigned int i = 0; i < hardBB.interpolant_.size(); ++i)
    {
        interpolant_[i] = new Interpolant(*hardBB.interpolant_[i]);
    }
}

HardBB::~HardBB()
{
    for(std::vector<Interpolant*>::iterator it = interpolant_.begin(); it != interpolant_.end(); ++it)
    {
        delete *it;
    }
}

bool HardBB::compare(const RealPhoton& photon) const
{
    const HardBB* bb = static_cast<const HardBB*>(&photon);

    if (interpolant_.size() != bb->interpolant_.size())
    {
        return false;
    }

    for (unsigned int i = 0; i < interpolant_.size(); ++i)
    {
        if (*interpolant_[i] != *bb->interpolant_[i])
            return false;
    }

    return RealPhoton::compare(photon);
}

// ------------------------------------------------------------------------- //
double HardBB::CalculateHardBB(double energy, double v)
{
    if (energy < 1.0e5 || v < 1.0e-7)
    {
        return 0;
    }

    double aux, sum, lov, loe;

    sum = 0;
    aux = 1;
    lov = log(v) / LOG10;
    loe = log(energy) / LOG10 - 3;

    for (unsigned int i = 0; i < interpolant_.size(); i++)
    {
        if (i > 0)
        {
            aux *= lov;
        }

        sum += aux * interpolant_[i]->InterpolateArray(loe);
    }
    return sum / v;
}

/******************************************************************************
*                                   SoftBB                                    *
******************************************************************************/

SoftBB::SoftBB()
{
}

SoftBB::SoftBB(const SoftBB& bb)
    : RealPhoton(bb)
{
}

SoftBB::~SoftBB()
{
}

double SoftBB::CalculateHardBB(double energy, double v)
{
    (void) energy;
    (void) v;

    return 0;
}


/******************************************************************************
*                                ShadowEffect                                *
******************************************************************************/

const std::string ShadowDuttaRenoSarcevicSeckel::name_ = "ShadowDRSS";
const std::string ShadowButkevichMikhailov::name_ = "ShadowButkevichMikhailov";

bool ShadowEffect::operator==(const ShadowEffect& shadow) const
{
    if (typeid(*this) != typeid(shadow))
        return false;
    else
        return true;
}

bool ShadowEffect::operator!=(const ShadowEffect& shadow) const
{
    return !(*this == shadow);
}

// ------------------------------------------------------------------------- //
// Dutta, Reno, Sarcevic, Seckel
// Phys Rev D 63 (2001), 094020
// eq. 3.10
// ------------------------------------------------------------------------- //
double ShadowDuttaRenoSarcevicSeckel::CalculateShadowEffect(const Components::Component& component, double x, double nu)
{
    (void) nu;

    if (component.GetNucCharge() == 1)
        return 1;

    if(x < 0.0014)
    {
        return pow(component.GetAtomicNum(), -0.1);
    }
    else if(x < 0.04)
    {
        return pow(component.GetAtomicNum(), 0.069 * log(x) / LOG10 + 0.097);
    }
    else
    {
        return 1;
    }
}

// ------------------------------------------------------------------------- //
size_t ShadowDuttaRenoSarcevicSeckel::GetHash() const
{
    size_t seed = 0;
    boost::hash_combine(seed, "ShadowDRSS");

    return seed;
}

// Butkevich, Mikheyev
// JETP 95 (2002), 11
// ------------------------------------------------------------------------- //
double ShadowButkevichMikhailov::CalculateShadowEffect(const Components::Component& component, double x, double nu)
{
    if(component.GetNucCharge()==1) return 1;

    double G;

    if(x > 0.3)
    {
        const double Mb = 0.437;
        const double la = 0.5;
        const double x2 = 0.278;

        double au = 1 / (1 - x);
        double ac = 1 / (1 - x2);
        // eq. 48
        double Aosc = (1 - la*x) * (au - ac
                - MPI / component.GetAverageNucleonWeight() * (au * au - ac * ac));
        // eq. 44
        G    = 1 - Mb * component.GetMN() * Aosc;
    }
    else
    {
        const double M1 = 0.129;
        const double M2 = 0.456;
        const double M3 = 0.553;

        double m1, m2, m3, x0, sgn, tmp;

        m1 = M1 * component.GetMN();
        m2 = M2 * component.GetMN();
        m3 = M3 * component.GetMN();
        nu *= 1.e-3;
        // eq. 53
        sgn = 112.2 * (0.609 * pow(nu, 0.0988) + 1.037 * pow(nu, -0.5944));

        // Bezrukav Bugaev shadow
        tmp = 0.00282*pow(component.GetAtomicNum(), 1. / 3) * sgn;
        G   = (3 /tmp) * (0.5 + ((1 + tmp) * exp(-tmp) - 1) / (tmp * tmp));

        // eq. 55
        G   = 0.75 * G + 0.25;
        x0  = pow(G / (1 + m2), 1 / m1);

        if(x>=x0)
        {
            // eq. 49
            G = pow(x, m1) * (1 + m2) * (1 - m3 * x);
        }
    }

    return G;
}

// ------------------------------------------------------------------------- //
size_t ShadowButkevichMikhailov::GetHash() const
{
    size_t seed = 0;
    boost::hash_combine(seed, "ShadowButkevichMikhailov");

    return seed;
}

/******************************************************************************
*                               Photonuclear                                  *
******************************************************************************/


// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

Photonuclear::Photonuclear(const ParticleDef& particle_def,
                           const Medium& medium,
                           const EnergyCutSettings& cuts,
                           double multiplier)
    : Parametrization(particle_def, medium, cuts, multiplier)
{
}

Photonuclear::Photonuclear(const Photonuclear& brems)
    : Parametrization(brems)
{
}

Photonuclear::~Photonuclear()
{
}

bool Photonuclear::compare(const Parametrization& parametrization) const
{
    return Parametrization::compare(parametrization);
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
Parametrization::IntegralLimits Photonuclear::GetIntegralLimits(double energy)
{
    double aux;

    IntegralLimits limits;

    limits.vMin = (MPI + MPI * MPI / (2 * components_[component_index_]->GetAverageNucleonWeight())) / energy;

    if(particle_def_.mass < MPI)
    {
        aux = particle_def_.mass/components_[component_index_]->GetAverageNucleonWeight();
        limits.vMax = 1 - components_[component_index_]->GetAverageNucleonWeight() * (1 + aux * aux) / (2 * energy);
    }
    else
    {
        limits.vMax = 1;
    }

    // vMax calculated above is always smaller than 1-m/E
    // in comparison, the following inequality arise
    // (M-m)^2 >= 0
    // limits.vMax = std::min(limits.vMax, 1 - particle_def_.mass/energy);

    if(limits.vMax < limits.vMin)
    {
        limits.vMax = limits.vMin;
    }

    limits.vUp = std::min(limits.vMax, cut_settings_.GetCut(energy));

    if(limits.vUp < limits.vMin)
    {
        limits.vUp = limits.vMin;
    }

    return limits;
}

