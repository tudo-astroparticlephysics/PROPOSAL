/*! \file   Medium.cxx
*   \brief  Source file for the medium routines.
*
*   For more details see the class documentation.
*
*   \date   2013.03.14
*   \author Jan-Hendrik Koehne
*/

// #include <cmath>
// #include <iomanip>

#include <boost/algorithm/string.hpp> // case insensitive string compare for configuration file

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/methods.h"

using namespace std;
using namespace PROPOSAL;

/******************************************************************************
*                                  OStream                                    *
******************************************************************************/

namespace PROPOSAL {

ostream& operator<<(ostream& os, Medium const& medium)
{

    os << "=================== Medium( " << &medium << " ) ===================" << endl;
    os << medium.name_ << endl;
    os << "\tnumber of components:\t\t\t\t" << medium.numComponents_ << endl;
    os << "\tmass density [g/cm3]:\t\t\t\t" << medium.massDensity_ << endl;
    os << "\tmolecule density [number/cm3]:\t\t\t" << medium.molDensity_ << endl;
    os << "\t<Z/A>:\t\t\t\t\t\t" << medium.ZA_ << endl;
    os << "\tsum of nucleons of all nuclei:\t\t\t" << medium.sumNucleons_ << endl;
    os << "\tionization potential [eV]:\t\t\t" << medium.I_ << endl;
    os << "\trefraction index:\t\t\t\t" << medium.r_ << endl;
    os << "\taverage all-component nucleon weight [MeV]:\t" << medium.MM_ << endl;
    os << "\tmultiplicative density correction factor:\t" << medium.rho_ << endl;
    os << "\tsum of charges of all nuclei:\t\t\t" << medium.sumCharge_ << endl;
    os << "\tradiation Length:\t\t" << medium.radiationLength_ << endl;

    for (std::vector<Components::Component*>::const_iterator iter = medium.components_.begin();
         iter != medium.components_.end();
         ++iter)
    {
        os << **iter << endl;
    }

    os << "=================================================================";
    return os;
}

} // namespace PROPOSAL

/******************************************************************************
*                 Temp. Functions for radiation length calc.                  *
******************************************************************************/

double fZ(unsigned int Z);
double Lrad(unsigned int Z);
double Lrad_dash(unsigned int Z);
double X0_inv(unsigned int Z, double M);

/******************************************************************************
*                                   Medium                                    *
******************************************************************************/

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
    , r0_(medium.r0_)
{
    // Deep copy of components
    components_.resize(numComponents_);
    for (unsigned int i = 0; i < components_.size(); ++i)
    {
        components_.at(i) = medium.components_.at(i)->clone();
    }
}

Medium::Medium(std::string name,
               double rho,
               double I,
               double C,
               double a,
               double m,
               double X0,
               double X1,
               double d0,
               double massDensity)
    : name_(name)
    , numComponents_(0)
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
    , r0_(0)
{
    if (rho > 0)
    {
        rho_ = rho;
    } else
    {
        rho_ = 1;
    }
}

// ------------------------------------------------------------------------- //
Medium::~Medium()
{
    for (std::vector<Components::Component*>::iterator iter = components_.begin();
         iter != components_.end();
         ++iter)
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
    swap(r0_, medium.r0_);

    components_.swap(medium.components_);
}

// ------------------------------------------------------------------------- //
Medium& Medium::operator=(const Medium& medium)
{
    if (this != &medium)
    {
        numComponents_ = medium.numComponents_;
        components_ = medium.components_;
        ZA_ = medium.ZA_;
        I_ = medium.I_;
        C_ = medium.C_;
        a_ = medium.a_;
        m_ = medium.m_;
        X0_ = medium.X0_;
        X1_ = medium.X1_;
        d0_ = medium.d0_;
        rho_ = medium.rho_;
        massDensity_ = medium.massDensity_;
        molDensity_ = medium.molDensity_;
        radiationLength_ = medium.radiationLength_;
        ecut_ = medium.ecut_;
        vcut_ = medium.vcut_;
        vCut_ = medium.vCut_;
        MM_ = medium.MM_;
        sumNucleons_ = medium.sumNucleons_;
        r0_ = medium.r0_;

        components_.clear();
        components_.resize(numComponents_);

        for (unsigned int i = 0; i < components_.size(); ++i)
        {
            components_.at(i) = medium.components_.at(i)->clone();
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
    else if (r0_ != medium.r0_)
        return false;
    else
    {
        bool Return = true;
        for (unsigned int i = 0; i < components_.size(); ++i)
        {
            if (*components_.at(i) != *medium.components_.at(i))
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

    for (std::vector<Components::Component*>::iterator iter = components_.begin();
         iter != components_.end();
         ++iter)
    {
        aux1 += (*iter)->GetAtomInMolecule() * (*iter)->GetNucCharge();
        aux2 += (*iter)->GetAtomInMolecule() * (*iter)->GetAtomicNum();
        aux3 += (*iter)->GetAtomInMolecule() * (*iter)->GetAtomicNum() *
                (*iter)->GetAverageNucleonWeight();
    }

    sumCharge_ = aux1;
    sumNucleons_ = aux2;
    ZA_ = aux1 / aux2;
    molDensity_ = massDensity_ * NA / aux2;
    MM_ = aux3 / aux2;
    r_ = 1.31; // only for ice - change if needed (sea water: 1.35)

    // Calculation of the radiation length
    aux1 = 0;
    aux2 = 0;

    for (std::vector<Components::Component*>::iterator iter = components_.begin();
         iter != components_.end();
         ++iter)
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

void Medium::SetComponents(std::vector<Components::Component*> components)
{
    components_ = components;
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

void Medium::SetR0(double r0)
{
    r0_ = r0;
}

/******************************************************************************
*                              Different Media                                *
******************************************************************************/

Water::Water(double rho)
    : Medium("Water",
             rho,
             75.0,    // I
             -3.5017, // C
             0.09116, // a
             3.4773,  // m
             0.2400,  // X0
             2.8004,  // X1
             0,       // d0
             1.000)   // massDensitiy
{
    components_.push_back(new Components::Hydrogen());
    components_.push_back(new Components::Oxygen());
    init();
}

Water& Water::operator=(const Water& medium)
{
    if (this != &medium) {
        const Water* med = dynamic_cast<const Water*>(&medium);
        if (!med) {
            log_warn("Cannot assign Water!");
            return *this;
        }

        Water tmp(*med);
        swap(tmp);
    }
    return *this;
}

Ice::Ice(double rho)
    : Medium("Ice",
             rho,
             75.0,    // I
             -3.5017, // C
             0.09116, // a
             3.4773,  // m
             0.2400,  // X0
             2.8004,  // X1
             0,       // d0
             0.917)   // massDensitiy
{
    components_.push_back(new Components::Hydrogen());
    components_.push_back(new Components::Oxygen());
    init();
}

Ice& Ice::operator=(const Ice& medium)
{
    if (this != &medium) {
        const Ice* med = dynamic_cast<const Ice*>(&medium);
        if (!med) {
            log_warn("Cannot assign Ice!");
            return *this;
        }

        Ice tmp(*med);
        swap(tmp);
    }
    return *this;
}

Salt::Salt(double rho)
    : Medium("Salt",
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
             2.323)   // Solid halite density
{
    components_.push_back(new Components::Natrium());
    components_.push_back(new Components::Chloride());
    init();
}

Salt& Salt::operator=(const Salt& medium)
{
    if (this != &medium) {
        const Salt* med = dynamic_cast<const Salt*>(&medium);
        if (!med) {
            log_warn("Cannot assign Salt!");
            return *this;
        }

        Salt tmp(*med);
        swap(tmp);
    }
    return *this;
}

StandardRock::StandardRock(double rho)
    : Medium("StandardRock",
             rho,
             136.4,   // I
             -3.7738, // C
             0.08301, // a
             3.4120,  // m
             0.0492,  // X0
             3.0549,  // X1
             0,       // d0
             2.650)   // massDensity
{
    components_.push_back(new Components::StandardRock());
    init();
}

StandardRock& StandardRock::operator=(const StandardRock& medium)
{
    if (this != &medium) {
        const StandardRock* med = dynamic_cast<const StandardRock*>(&medium);
        if (!med) {
            log_warn("Cannot assign StandardRock!");
            return *this;
        }

        StandardRock tmp(*med);
        swap(tmp);
    }
    return *this;
}

FrejusRock::FrejusRock(double rho)
    : Medium("FrejusRock",
             rho,
             149.0,  // I
             -5.053, // C
             0.078,  // a
             3.645,  // m
             0.288,  // X0
             3.196,  // X1
             0,      // d0
             2.740)  // massDensity
{
    components_.push_back(new Components::FrejusRock());
    init();
}

FrejusRock& FrejusRock::operator=(const FrejusRock& medium)
{
    if (this != &medium) {
        const FrejusRock* med = dynamic_cast<const FrejusRock*>(&medium);
        if (!med) {
            log_warn("Cannot assign FrejusRock!");
            return *this;
        }

        FrejusRock tmp(*med);
        swap(tmp);
    }
    return *this;
}

Iron::Iron(double rho)
    : Medium("Iron",
             rho,
             286.0,   // I
             -4.2911, // C
             0.14680, // a
             2.9632,  // m
             -0.0012, // X0
             3.1531,  // X1
             0.12,    // d0
             7.874)   // massDensity
{
    components_.push_back(new Components::Iron());
    init();
}

Iron& Iron::operator=(const Iron& medium)
{
    if (this != &medium) {
        const Iron* med = dynamic_cast<const Iron*>(&medium);
        if (!med) {
            log_warn("Cannot assign Iron!");
            return *this;
        }

        Iron tmp(*med);
        swap(tmp);
    }
    return *this;
}

Hydrogen::Hydrogen(double rho)
    : Medium("Hydrogen",
             rho,
             21.8,    // I
             -3.0977, // C
             0.13483, // a
             5.6249,  // m
             0.4400,  // X0
             1.8856,  // X1
             0,       // d0
             0.07080) // massDensity
{
    components_.push_back(new Components::Hydrogen());
    init();
}

Hydrogen& Hydrogen::operator=(const Hydrogen& medium)
{
    if (this != &medium) {
        const Hydrogen* med = dynamic_cast<const Hydrogen*>(&medium);
        if (!med) {
            log_warn("Cannot assign Hydrogen!");
            return *this;
        }

        Hydrogen tmp(*med);
        swap(tmp);
    }
    return *this;
}

Lead::Lead(double rho)
    : Medium("Lead",
             rho,
             823.0,   // I
             -6.2018, // C
             0.09359, // a
             3.1608,  // m
             0.3776,  // X0
             3.8073,  // X1
             0.14,    // d0
             11.350)  // massDensity
{
    components_.push_back(new Components::Lead());
    init();
}

Lead& Lead::operator=(const Lead& medium)
{
    if (this != &medium) {
        const Lead* med = dynamic_cast<const Lead*>(&medium);
        if (!med) {
            log_warn("Cannot assign Lead!");
            return *this;
        }

        Lead tmp(*med);
        swap(tmp);
    }
    return *this;
}

Copper::Copper(double rho)
    : Medium("Copper",
             rho,
             322.0,   // I
             -4.4190, // C
             0.14339, // a
             2.9044,  // m
             -0.0254, // X0
             3.2792,  // X1
             0.08,    // d0
             8.960)   // massDensity
{
    components_.push_back(new Components::Copper());
    init();
}

Copper& Copper::operator=(const Copper& medium)
{
    if (this != &medium) {
        const Copper* med = dynamic_cast<const Copper*>(&medium);
        if (!med) {
            log_warn("Cannot assign Copper!");
            return *this;
        }

        Copper tmp(*med);
        swap(tmp);
    }
    return *this;
}

Uranium::Uranium(double rho)
    : Medium("Uranium",
             rho,
             890.0,   // I
             -5.8694, // C
             0.19677, // a
             2.8171,  // m
             0.2260,  // X0
             3.3721,  // X1
             0.14,    // d0
             18.950)  // massDensity
{
    components_.push_back(new Components::Uranium());
    init();
}

Uranium& Uranium::operator=(const Uranium& medium)
{
if (this != &medium) {
    const Uranium* med = dynamic_cast<const Uranium*>(&medium);
    if (!med) {
        log_warn("Cannot assign Uranium!");
        return *this;
    }

    Uranium tmp(*med);
    swap(tmp);
}
return *this;
}

Air::Air(double rho)
    : Medium("Air",
             rho,
             85.7,     // I
             -10.5961, // C
             0.10914,  // a
             3.3994,   // m
             1.7418,   // X0
             4.2759,   // X1
             0,        // d0
             1.205e-3) // dry, 1 atm massDensity
{
    const double fr1 = 2 * 78.1;
    const double fr2 = 2 * 21.0;
    const double fr3 = 0.9;
    const double fra = fr1 + fr2 + fr3;

    components_.push_back(new Components::Nitrogen(fr1 / fra));
    components_.push_back(new Components::Oxygen(fr2 / fra));
    components_.push_back(new Components::Arsenic(fr3 / fra));
    init();
}

Air& Air::operator=(const Air& medium)
{
    if (this != &medium) {
        const Air* med = dynamic_cast<const Air*>(&medium);
        if (!med) {
            log_warn("Cannot assign Air!");
            return *this;
        }

        Air tmp(*med);
        swap(tmp);
    }
    return *this;
}

Paraffin::Paraffin(double rho)
    : Medium("Paraffin",
             rho,
             55.9,    // I
             -2.9551, // C
             0.1209,  // a
             3.4288,  // m
             0.1289,  // X0
             2.5084,  // X1
             0,       // d0
             0.93)    // massDensity
{
    components_.push_back(new Components::Carbon(25.0));
    components_.push_back(new Components::Hydrogen(52.0));
    init();
}

Paraffin& Paraffin::operator=(const Paraffin& medium)
{
    if (this != &medium) {
        const Paraffin* med = dynamic_cast<const Paraffin*>(&medium);
        if (!med) {
            log_warn("Cannot assign Paraffin!");
            return *this;
        }

        Paraffin tmp(*med);
        swap(tmp);
    }
    return *this;
}

AntaresWater::AntaresWater(double rho)
    : Medium("AntaresWater",
             rho,
             75.0,    // I
             -3.5017, // C
             0.09116, // a
             3.4773,  // m
             0.2400,  // X0
             2.8004,  // X1
             0,       // d0
             1.03975) // massDensity
{
    //  Chemical composition of the seawater
    //  according to
    //  A.Okada, Astropart. Phys. 2 (1994) 393
    //  and references therein
    //  corrected for Mediterranean Sea, ANTARES place
    //  according to salinity  38.44+-0.02 g/kg,
    //  as cited in J.Brunner, ANTARES-Site/2000-001
    //  instead of 35.0 g/kg as cited in A.Okada, ...
    //  (so, n[2-7] have been just multiplied by 1.098)

    components_.push_back(new Components::Hydrogen(2.0));
    components_.push_back(new Components::Oxygen(1.00884));
    components_.push_back(new Components::Natrium(0.00943));
    components_.push_back(new Components::Potassium(0.000209));
    components_.push_back(new Components::Magnesium(0.001087));
    components_.push_back(new Components::Calcium(0.000209));
    components_.push_back(new Components::Chloride(0.01106));
    components_.push_back(new Components::Sulfur(0.00582));
    init();

    // J.Brunner, ANTARES-Site/2000-001, the mean value
    // for sea water density at the ANTARES place between
    // sea bed D = 2400 m (1.0404 g/cm^3) and middle of
    // detector D = 2126 m (1.0391 g/cm^3).
    // Ro = 1.0341
    // for sea water density at the ANTARES place between
    // surface D = 0 m (1.0291 g/cm^3) and middle of
    // detector D = 2126 m (1.0391 g/cm^3)
}

AntaresWater& AntaresWater::operator=(const AntaresWater& medium)
{
    if (this != &medium) {
        const AntaresWater* med = dynamic_cast<const AntaresWater*>(&medium);
        if (!med) {
            log_warn("Cannot assign AntaresWater!");
            return *this;
        }

        AntaresWater tmp(*med);
        swap(tmp);
    }
    return *this;
}

/******************************************************************************
*                               Medium Factory                                *
******************************************************************************/

MediumFactory::MediumFactory()
{
    // Register all media in lower case!

    Register("water", &Water::create);
    Register("ice", &Ice::create);
    Register("salt", &Salt::create);
    Register("standardrock", &StandardRock::create);
    Register("frejusrock", &FrejusRock::create);
    Register("iron", &Iron::create);
    Register("hydrogen", &Hydrogen::create);
    Register("lead", &Lead::create);
    Register("copper", &Copper::create);
    Register("uranium", &Uranium::create);
    Register("air", &Air::create);
    Register("paraffin", &AntaresWater::create);
}

void MediumFactory::Register(const std::string& name, Medium* (*create)(void))
{
    medium_map[name] = create;
}

Medium* MediumFactory::CreateMedium(const std::string& name)
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    std::map<std::string, Medium* (*)(void)>::iterator it = medium_map.find(name_lower);

    if (it != medium_map.end())
    {
        return it->second();
    } else
    {
        return NULL;
    }
}

/******************************************************************************
*                              Helper Funcitons                               *
******************************************************************************/

double fZ(unsigned int Z)
{
    double a_sq = ALPHA * ALPHA * Z * Z;

    return a_sq * (1. / (1. + a_sq) + 0.20206 - 0.0369 * a_sq + 0.0083 * a_sq * a_sq -
                   0.002 * a_sq * a_sq * a_sq);
}

double Lrad(unsigned int Z)
{
    if (Z > 4)
        return log(184.15 * pow(Z, -1. / 3.)); // Elements Z>4

    if (Z == 1)
        return 5.31; // Hydrogen
    if (Z == 2)
        return 4.79; // Helium
    if (Z == 3)
        return 4.74; // Lithium
    if (Z == 4)
        return 4.71; // Beryllium
    return 0;
}

double Lrad_dash(unsigned int Z)
{
    if (Z > 4)
        return log(1194. * pow(Z, -2. / 3.)); // Elements Z>4

    if (Z == 1)
        return 6.144; // Hydrogen
    if (Z == 2)
        return 5.621; // Helium
    if (Z == 3)
        return 5.805; // Lithium
    if (Z == 4)
        return 5.924; // Beryllium
    return 0;
}

double X0_inv(unsigned int Z, double M)
{
    return 4. * ALPHA * RE * RE * NA / M * (Z * Z * (Lrad(Z) - fZ(Z)) + Z * Lrad_dash(Z));
}
