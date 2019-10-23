/*! \file   Medium.cxx
 *   \brief  Source file for the medium routines.
 *
 *   For more details see the class documentation.
 *
 *   \date   2013.03.14
 *   \author Jan-Hendrik Koehne
 */

#include <cmath>
#include <sstream>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

/******************************************************************************
 *                                  OStream                                    *
 ******************************************************************************/

namespace PROPOSAL {

std::ostream& operator<<(std::ostream& os, Medium const& medium)
{

    std::stringstream ss;
    ss << " Medium (" << &medium << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    os << medium.name_ << std::endl;
    os << "number of components:\t\t\t\t" << medium.numComponents_ << std::endl;
    os << "mass density [g/cm3]:\t\t\t\t" << medium.massDensity_ << std::endl;
    os << "molecule density [number/cm3]:\t\t\t" << medium.molDensity_ << std::endl;
    os << "<Z/A>:\t\t\t\t\t\t" << medium.ZA_ << std::endl;
    os << "sum of nucleons of all nuclei:\t\t\t" << medium.sumNucleons_ << std::endl;
    os << "ionization potential [eV]:\t\t\t" << medium.I_ << std::endl;
    os << "refraction index:\t\t\t\t" << medium.r_ << std::endl;
    os << "average all-component nucleon weight [MeV]:\t" << medium.MM_ << std::endl;
    os << "multiplicative density correction factor:\t" << medium.rho_ << std::endl;
    os << "sum of charges of all nuclei:\t\t\t" << medium.sumCharge_ << std::endl;
    os << "radiation Length:\t\t" << medium.radiationLength_ << std::endl;

    for (std::vector<Components::Component*>::const_iterator iter = medium.components_.begin();
         iter != medium.components_.end();
         ++iter)
    {
        os << **iter << std::endl;
    }

    os << Helper::Centered(60, "");
    return os;
}

} // namespace PROPOSAL

/******************************************************************************
 *                                   Medium                                    *
 ******************************************************************************/

Medium::Medium(std::string name,
               double rho,
               double I,
               double C,
               double a,
               double m,
               double X0,
               double X1,
               double d0,
               double massDensity,
               std::vector<std::shared_ptr<Components::Component>> components)
    : name_(name)
    , numComponents_(components.size())
    , sumCharge_(0)
    , ZA_(0)
    , I_(I)
    , C_(C)
    , a_(a)
    , m_(m)
    , X0_(X0)
    , X1_(X1)
    , d0_(d0)
    , r_(0)
    , massDensity_(massDensity)
    , molDensity_(0)
    , radiationLength_(0)
    , ecut_(0)
    , vcut_(0)
    , vCut_(0)
    , MM_(0)
    , sumNucleons_(0)
{
    if (rho > 0)
    {
        rho_ = rho;
    } else
    {
        rho_ = 1;
    }

    // components_.reserve(numComponents_);
    // for (int i = 0; i < numComponents_; ++i)
    for (auto component: components)
    {
        components_.push_back(component->clone());
    }

    init();
}

Medium::Medium(const Medium& medium)
    : name_(medium.name_)
    , numComponents_(medium.numComponents_)
    , sumCharge_(medium.sumCharge_)
    , ZA_(medium.ZA_)
    , I_(medium.I_)
    , C_(medium.C_)
    , a_(medium.a_)
    , m_(medium.m_)
    , X0_(medium.X0_)
    , X1_(medium.X1_)
    , d0_(medium.d0_)
    , r_(medium.r_)
    , rho_(medium.rho_)
    , massDensity_(medium.massDensity_)
    , molDensity_(medium.molDensity_)
    , radiationLength_(medium.radiationLength_)
    , ecut_(medium.ecut_)
    , vcut_(medium.vcut_)
    , vCut_(medium.vCut_)
    , MM_(medium.MM_)
    , sumNucleons_(medium.sumNucleons_)
{
    // Deep copy of components
    components_.resize(numComponents_);
    for (unsigned int i = 0; i < components_.size(); ++i)
    {
        components_[i] = medium.components_[i]->clone();
    }
}

// ------------------------------------------------------------------------- //
Medium::~Medium()
{
    for (std::vector<Components::Component*>::iterator iter = components_.begin(); iter != components_.end(); ++iter)
    {
        delete (*iter);
    }

    components_.clear();
}

// ------------------------------------------------------------------------- //
// Operators & swap
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
void Medium::swap(Medium& medium)
{
    using std::swap;

    swap(name_, medium.name_);
    swap(numComponents_, medium.numComponents_);
    swap(sumCharge_, medium.sumCharge_);
    swap(ZA_, medium.ZA_);
    swap(I_, medium.I_);
    swap(C_, medium.C_);
    swap(a_, medium.a_);
    swap(m_, medium.m_);
    swap(X0_, medium.X0_);
    swap(X1_, medium.X1_);
    swap(d0_, medium.d0_);
    swap(r_, medium.r_);
    swap(rho_, medium.rho_);
    swap(massDensity_, medium.massDensity_);
    swap(molDensity_, medium.molDensity_);
    swap(radiationLength_, medium.radiationLength_);
    swap(ecut_, medium.ecut_);
    swap(vcut_, medium.vcut_);
    swap(vCut_, medium.vCut_);
    swap(MM_, medium.MM_);
    swap(sumNucleons_, medium.sumNucleons_);

    components_.swap(medium.components_);
}

// ------------------------------------------------------------------------- //
Medium& Medium::operator=(const Medium& medium)
{
    if (this != &medium)
    {
        numComponents_   = medium.numComponents_;
        components_      = medium.components_;
        ZA_              = medium.ZA_;
        I_               = medium.I_;
        C_               = medium.C_;
        a_               = medium.a_;
        m_               = medium.m_;
        X0_              = medium.X0_;
        X1_              = medium.X1_;
        d0_              = medium.d0_;
        rho_             = medium.rho_;
        massDensity_     = medium.massDensity_;
        molDensity_      = medium.molDensity_;
        radiationLength_ = medium.radiationLength_;
        ecut_            = medium.ecut_;
        vcut_            = medium.vcut_;
        vCut_            = medium.vCut_;
        MM_              = medium.MM_;
        sumNucleons_     = medium.sumNucleons_;

        components_.clear();
        components_.resize(numComponents_);

        for (unsigned int i = 0; i < components_.size(); ++i)
        {
            components_[i] = medium.components_[i]->clone();
        }
    }

    return *this;
}

// ------------------------------------------------------------------------- //
bool Medium::operator==(const Medium& medium) const
{
    if (name_ != medium.name_)
        return false;
    else if (numComponents_ != medium.numComponents_)
        return false;
    else if (sumCharge_ != medium.sumCharge_)
        return false;
    else if (ZA_ != medium.ZA_)
        return false;
    else if (I_ != medium.I_)
        return false;
    else if (C_ != medium.C_)
        return false;
    else if (a_ != medium.a_)
        return false;
    else if (m_ != medium.m_)
        return false;
    else if (X0_ != medium.X0_)
        return false;
    else if (X1_ != medium.X1_)
        return false;
    else if (d0_ != medium.d0_)
        return false;
    else if (r_ != medium.r_)
        return false;
    else if (rho_ != medium.rho_)
        return false;
    else if (massDensity_ != medium.massDensity_)
        return false;
    else if (molDensity_ != medium.molDensity_)
        return false;
    else if (radiationLength_ != medium.radiationLength_)
        return false;
    else if (ecut_ != medium.ecut_)
        return false;
    else if (vcut_ != medium.vcut_)
        return false;
    else if (vCut_ != medium.vCut_)
        return false;
    else if (MM_ != medium.MM_)
        return false;
    else if (sumNucleons_ != medium.sumNucleons_)
        return false;
    else
    {
        bool Return = true;
        for (unsigned int i = 0; i < components_.size(); ++i)
        {
            if (*components_[i] != *medium.components_[i])
            {
                Return = false;
            }
        }

        return Return;
    }
}

// ------------------------------------------------------------------------- //
bool Medium::operator!=(const Medium& medium) const
{
    return !(*this == medium);
}

// ------------------------------------------------------------------------- //
// Methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
void Medium::init()
{
    // Init Radiation length and further members
    double aux1 = 0;
    double aux2 = 0;
    double aux3 = 0;

    numComponents_ = components_.size();
    massDensity_ *= rho_;

    for (std::vector<Components::Component*>::iterator iter = components_.begin(); iter != components_.end(); ++iter)
    {
        aux1 += (*iter)->GetAtomInMolecule() * (*iter)->GetNucCharge();
        aux2 += (*iter)->GetAtomInMolecule() * (*iter)->GetAtomicNum();
        aux3 += (*iter)->GetAtomInMolecule() * (*iter)->GetAtomicNum() * (*iter)->GetAverageNucleonWeight();
    }

    sumCharge_   = aux1;
    sumNucleons_ = aux2;
    ZA_          = aux1 / aux2;
    molDensity_  = massDensity_ * NA / aux2;
    MM_          = aux3 / aux2;

    // TODO: this is never used; is that art or deletable
    r_ = 1.31; // only for ice - change if needed (sea water: 1.35)

    // TODO: Compare to Bremsstrahlung::CalculateScatteringX0; just one (this or the Brems-thing) is needed
    // Calculation of the radiation length
    aux1 = 0;
    aux2 = 0;

    for (std::vector<Components::Component*>::iterator iter = components_.begin(); iter != components_.end(); ++iter)
    {
        aux1 += X0_inv((*iter)->GetNucCharge(), (*iter)->GetAtomicNum()) *
                ((*iter)->GetAtomInMolecule() * (*iter)->GetAtomicNum());
        aux2 += (*iter)->GetAtomInMolecule() * (*iter)->GetAtomicNum();
    }

    radiationLength_ = aux2 / aux1;
    radiationLength_ /= massDensity_;
}

// ------------------------------------------------------------------------- //
// Setter
// ------------------------------------------------------------------------- //

void Medium::SetComponents(std::vector<std::shared_ptr<Components::Component>> components)
{
    for (auto component: components_)
    {
        delete component;
    }

    components_.clear();

    for (auto component: components)
    {
        components_.push_back(component->clone());
    }

    init(); // Init further member according to these components
}

void Medium::SetZA(double ZA)
{
    ZA_ = ZA;
}

void Medium::SetI(double I)
{
    I_ = I;
}

void Medium::SetC(double C)
{
    C_ = C;
}

void Medium::SetA(double a)
{
    a_ = a;
}

void Medium::SetM(double m)
{
    m_ = m;
}

void Medium::SetX0(double X0)
{
    X0_ = X0;
}

void Medium::SetX1(double X1)
{
    X1_ = X1;
}

void Medium::SetD0(double d0)
{
    d0_ = d0;
}

void Medium::SetR(double r)
{
    r_ = r;
}

void Medium::SetRho(double rho)
{
    rho_ = rho;
}

void Medium::SetMassDensity(double massDensity)
{
    massDensity_ = massDensity;
}

void Medium::SetMolDensity(double molDensity)
{
    molDensity_ = molDensity;
}

void Medium::SetMM(double MM)
{
    MM_ = MM;
}

void Medium::SetSumNucleons(double sumNucleons)
{
    sumNucleons_ = sumNucleons;
}

/******************************************************************************
 *                              Different Media                                *
 ******************************************************************************/

Water::Water(double rho)
    : Medium("water",
             rho,
             75.0,    // I
             -3.5017, // C
             0.09116, // a
             3.4773,  // m
             0.2400,  // X0
             2.8004,  // X1
             0,       // d0
             1.000,   // massDensitiy
             {
                 std::shared_ptr<Components::Component>(new Components::Hydrogen(2)),
                 std::shared_ptr<Components::Component>(new Components::Oxygen())
             })
{
}

Ice::Ice(double rho)
    : Medium("ice",
             rho,
             75.0,    // I
             -3.5017, // C
             0.09116, // a
             3.4773,  // m
             0.2400,  // X0
             2.8004,  // X1
             0,       // d0
             0.917,   // massDensitiy
             {
                 std::shared_ptr<Components::Component>(new Components::Hydrogen(2)),
                 std::shared_ptr<Components::Component>(new Components::Oxygen())
             })
{
}

Salt::Salt(double rho)
    : Medium("salt",
             rho,
             // Calculated by ESTAR detabase
             // (it could be 185 eV by the method of reference below)
             175.3,   // I
             -4.5041, // C
             0.1632,  // a
             3,       // m
             0.2,     // X0
             3.0,     // X1
             0,       // d0
             2.323,   // Solid halite density
             {
                 std::shared_ptr<Components::Component>(new Components::Sodium()),
                 std::shared_ptr<Components::Component>(new Components::Chlorine())
             })
{
}

CalciumCarbonate::CalciumCarbonate(double rho)
    : Medium("calciumcarbonate",
             rho,
             136.4,   // I
             -3.7738, // C
             0.08301, // a
             3.4120,  // m
             0.0492,  // X0
             3.0549,  // X1
             0,       // d0
             2.650,   // massDensity
             {
                 std::shared_ptr<Components::Component>(new Components::Calcium()),
                 std::shared_ptr<Components::Component>(new Components::Carbon()),
                 std::shared_ptr<Components::Component>(new Components::Oxygen(3))
             })
{
}

StandardRock::StandardRock(double rho)
    : Medium("standardrock",
             rho,
             136.4,   // I
             -3.7738, // C
             0.08301, // a
             3.4120,  // m
             0.0492,  // X0
             3.0549,  // X1
             0,       // d0
             2.650,   // massDensity
             {std::shared_ptr<Components::Component>(new Components::StandardRock())})
{
}

FrejusRock::FrejusRock(double rho)
    : Medium("frejusrock",
             rho,
             149.0,  // I
             -5.053, // C
             0.078,  // a
             3.645,  // m
             0.288,  // X0
             3.196,  // X1
             0,      // d0
             2.740,  // massDensity
             {std::shared_ptr<Components::Component>(new Components::FrejusRock())})
{
}

Iron::Iron(double rho)
    : Medium("iron",
             rho,
             286.0,   // I
             -4.2911, // C
             0.14680, // a
             2.9632,  // m
             -0.0012, // X0
             3.1531,  // X1
             0.12,    // d0
             7.874,   // massDensity
             {std::shared_ptr<Components::Component>(new Components::Iron())})
{
}

Hydrogen::Hydrogen(double rho)
    : Medium("hydrogen",
             rho,
             21.8,    // I
             -3.0977, // C
             0.13483, // a
             5.6249,  // m
             0.4400,  // X0
             1.8856,  // X1
             0,       // d0
             0.07080, // massDensity
             {std::shared_ptr<Components::Component>(new Components::Hydrogen())})
{
}

Lead::Lead(double rho)
    : Medium("lead",
             rho,
             823.0,   // I
             -6.2018, // C
             0.09359, // a
             3.1608,  // m
             0.3776,  // X0
             3.8073,  // X1
             0.14,    // d0
             11.350,  // massDensity
             {std::shared_ptr<Components::Component>(new Components::Lead())})
{
}

Copper::Copper(double rho)
    : Medium("copper",
             rho,
             322.0,   // I
             -4.4190, // C
             0.14339, // a
             2.9044,  // m
             -0.0254, // X0
             3.2792,  // X1
             0.08,    // d0
             8.960,   // massDensity
             {std::shared_ptr<Components::Component>(new Components::Copper())})
{
}

Uranium::Uranium(double rho)
    : Medium("uranium",
             rho,
             890.0,   // I
             -5.8694, // C
             0.19677, // a
             2.8171,  // m
             0.2260,  // X0
             3.3721,  // X1
             0.14,    // d0
             18.950,  // massDensity
             {std::shared_ptr<Components::Component>(new Components::Uranium())})
{
}

const double Air::fraction_N   = 2 * 78.1;
const double Air::fraction_O   = 2 * 21.0;
const double Air::fraction_Ar  = 0.9;
const double Air::fraction_sum = Air::fraction_N + Air::fraction_O + Air::fraction_Ar;

Air::Air(double rho)
    : Medium("air",
             rho,
             85.7,     // I
             -10.5961, // C
             0.10914,  // a
             3.3994,   // m
             1.7418,   // X0
             4.2759,   // X1
             0,        // d0
             1.205e-3, // dry, 1 atm massDensity
             {
                 std::shared_ptr<Components::Component>(new Components::Nitrogen(fraction_N / fraction_sum)),
                 std::shared_ptr<Components::Component>(new Components::Oxygen(fraction_O / fraction_sum)),
                 std::shared_ptr<Components::Component>(new Components::Argon(fraction_Ar / fraction_sum))
             })
{
}

Paraffin::Paraffin(double rho)
    : Medium(
          "paraffin",
          rho,
          55.9,    // I
          -2.9551, // C
          0.1209,  // a
          3.4288,  // m
          0.1289,  // X0
          2.5084,  // X1
          0,       // d0
          0.93,    // massDensity
         {
             std::shared_ptr<Components::Component>(new Components::Carbon(25.0)),
             std::shared_ptr<Components::Component>(new Components::Oxygen(52.0))
         })
{
}

AntaresWater::AntaresWater(double rho)
    : Medium("antareswater",
             rho,
             75.0,    // I
             -3.5017, // C
             0.09116, // a
             3.4773,  // m
             0.2400,  // X0
             2.8004,  // X1
             0,       // d0
             1.03975, // massDensity
             {
                 std::shared_ptr<Components::Component>(new Components::Hydrogen(2.0)),
                 std::shared_ptr<Components::Component>(new Components::Oxygen(1.00884)),
                 std::shared_ptr<Components::Component>(new Components::Sodium(0.00943)),
                 std::shared_ptr<Components::Component>(new Components::Potassium(0.000209)),
                 std::shared_ptr<Components::Component>(new Components::Magnesium(0.001087)),
                 std::shared_ptr<Components::Component>(new Components::Calcium(0.000209)),
                 std::shared_ptr<Components::Component>(new Components::Chlorine(0.01106)),
                 std::shared_ptr<Components::Component>(new Components::Sulfur(0.00582))
             })
{
    /*
    * initialize ANTARES water
    * Sea water (Mediterranean Sea, ANTARES place)
    * ==========================================================================
    * WATER DENSITY CHANGES WITH THE DEPTH FROM 1.0291 g/cm^3 AT SURFACE
    * UP TO 1.0404 g/cm^3 AT THE SEA BED
    * (ANTARES-Site/2000-001 and references therein)
    *
    * The error which is caused by this simplified approach (average value for
    * density) does not exceed 0.5% (much less, in fact) that is comparable with
    *  an error which comes from uncertainties with the muon cross-sections.
    *==========================================================================
    */
    //  Chemical composition of the seawater
    //  according to
    //  A.Okada, Astropart. Phys. 2 (1994) 393
    //  and references therein
    //  corrected for Mediterranean Sea, ANTARES place
    //  according to salinity  38.44+-0.02 g/kg,
    //  as cited in J.Brunner, ANTARES-Site/2000-001
    //  instead of 35.0 g/kg as cited in A.Okada, ...
    //  (so, n[2-7] have been just multiplied by 1.098)

    // J.Brunner, ANTARES-Site/2000-001, the mean value
    // for sea water density at the ANTARES place between
    // sea bed D = 2400 m (1.0404 g/cm^3) and middle of
    // detector D = 2126 m (1.0391 g/cm^3).
    // Ro = 1.0341
    // for sea water density at the ANTARES place between
    // surface D = 0 m (1.0291 g/cm^3) and middle of
    // detector D = 2126 m (1.0391 g/cm^3)

}

// chemical composition is normalized with respect to hydrogen: calculate a
// total atomic weight that yields two hydrogen atoms; use this weight to
// calculate the numbers of atoms for all other elements
CascadiaBasinWater::CascadiaBasinWater(double rho)
    : Medium(
        "cascadiabasinwater",
        rho,
        75.,       // I
        -3.5017,   // C
        0.09116,   // a
        3.4773,    // m
        0.2400,    // X0
        2.8004,    // X1
        0.,        // d0
        1.0400322, // massDensity
        {
            std::make_shared<Components::Component>(Components::Hydrogen(2.0)),
            std::make_shared<Components::Component>(Components::Oxygen(1.0021)),
            std::make_shared<Components::Component>(Components::Sodium(8.7174e-3)),
            std::make_shared<Components::Component>(Components::Magnesium(9.8180e-4)),
            std::make_shared<Components::Component>(Components::Calcium(1.9113e-4)),
            std::make_shared<Components::Component>(Components::Potassium(1.8975e-4)),
            std::make_shared<Components::Component>(Components::Chlorine(1.0147e-2)),
            std::make_shared<Components::Component>(Components::Sulfur(5.2487e-4))
        })
{}

/******************************************************************************
 *                        private Helper Funcitons                             *
 ******************************************************************************/

double Medium::X0_inv(unsigned int Z, double M)
{
    double a_sq      = 0.;
    double fZ        = 0.;
    double Lrad      = 0.;
    double Lrad_dash = 0.;

    a_sq = ALPHA * ALPHA * Z * Z;
    fZ   = a_sq * (1. / (1. + a_sq) + 0.20206 - 0.0369 * a_sq + 0.0083 * a_sq * a_sq - 0.002 * a_sq * a_sq * a_sq);

    // Get radiation logarithm from Tsai (Rev.Mod.Phys.46)
    if (Z > 4) // Elements Z>4
    {
        // Thomas-Fermi-Moliere model (Tsai eq. B21)
        Lrad      = std::log(184.15 * std::pow(Z, -1. / 3.));
        Lrad_dash = std::log(1194. * std::pow(Z, -2. / 3.));
    }
    // Tsai table B2
    else if (Z == 1) // Hydrogen
    {
        Lrad      = 5.31;
        Lrad_dash = 6.144;
    } else if (Z == 2) // Helium
    {
        Lrad      = 4.79;
        Lrad_dash = 5.621;
    } else if (Z == 3) // Lithium
    {
        Lrad      = 4.74;
        Lrad_dash = 5.805;
    } else if (Z == 4) // Beryllium
    {
        Lrad      = 4.71;
        Lrad_dash = 5.924;
    }
    // Bremsstrahlung complete screening case
    return 4. * ALPHA * RE * RE * NA / M * (Z * Z * (Lrad - fZ) + Z * Lrad_dash);
}
