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
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

std::unique_ptr<std::map<size_t, Medium>> Medium::medium_map = nullptr;

/******************************************************************************
 *                                  OStream                                    *
 ******************************************************************************/

namespace PROPOSAL {

std::ostream& operator<<(std::ostream& os, Medium const& medium) {
    std::stringstream ss;
    ss << " Medium (" << &medium << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    os << medium.name_ << std::endl;
    os << "number of components:\t\t\t\t" << medium.numComponents_ << std::endl;
    os << "mass density [g/cm3]:\t\t\t\t" << medium.massDensity_ << std::endl;
    os << "molecule density [number/cm3]:\t\t\t" << medium.molDensity_
       << std::endl;
    os << "<Z/A>:\t\t\t\t\t\t" << medium.ZA_ << std::endl;
    os << "sum of nucleons of all nuclei:\t\t\t" << medium.sumNucleons_
       << std::endl;
    os << "ionization potential [eV]:\t\t\t" << medium.I_ << std::endl;
    os << "average all-component nucleon weight [MeV]:\t" << medium.MM_
       << std::endl;
    os << "sum of charges of all nuclei:\t\t\t" << medium.sumCharge_
       << std::endl;
    os << "radiation Length:\t\t" << medium.radiationLength_ << std::endl;

    for (std::vector<Component>::const_iterator iter =
             medium.components_.begin();
         iter != medium.components_.end(); ++iter) {
        os << *iter << std::endl;
    }

    os << Helper::Centered(60, "");
    return os;
}

}  // namespace PROPOSAL

/******************************************************************************
 *                                   Medium                                    *
 ******************************************************************************/

Medium::Medium(std::string name,
               double I,
               double C,
               double a,
               double m,
               double X0,
               double X1,
               double d0,
               double massDensity,
               std::vector<Component> components)
    : name_(name),
      numComponents_(components.size()),
      sumCharge_(0),
      ZA_(0),
      I_(I),
      C_(C),
      a_(a),
      m_(m),
      X0_(X0),
      X1_(X1),
      d0_(d0),
      massDensity_(massDensity),
      molDensity_(0),
      radiationLength_(0),
      MM_(0),
      sumNucleons_(0),
      components_(components) {
    init();

    if (!medium_map)
        medium_map = std::make_unique<std::map<size_t, Medium>>();

    auto hash = GetHash();
    if (medium_map->find(hash) == medium_map->end())
        medium_map->insert({hash, Medium(*this)});

}

// ------------------------------------------------------------------------- //
bool Medium::operator==(const Medium& medium) const {
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
    else if (massDensity_ != medium.massDensity_)
        return false;
    else if (molDensity_ != medium.molDensity_)
        return false;
    else if (radiationLength_ != medium.radiationLength_)
        return false;
    else if (MM_ != medium.MM_)
        return false;
    else if (sumNucleons_ != medium.sumNucleons_)
        return false;
    else {
        bool Return = true;
        for (unsigned int i = 0; i < components_.size(); ++i) {
            if (!(components_[i] == medium.components_[i])) {
                Return = false;
            }
        }

        return Return;
    }
}

bool Medium::operator!=(const Medium& medium) const
{
    return !(*this == medium);
}

// ------------------------------------------------------------------------- //
//                                   Methods
// ------------------------------------------------------------------------- //

void Medium::init() {
    // Init Radiation length and further members
    double aux1 = 0;
    double aux2 = 0;
    double aux3 = 0;

    numComponents_ = components_.size();

    for (std::vector<Component>::iterator iter =
             components_.begin();
         iter != components_.end(); ++iter) {
        aux1 += (iter)->GetAtomInMolecule() * (iter)->GetNucCharge();
        aux2 += (iter)->GetAtomInMolecule() * (iter)->GetAtomicNum();
        aux3 += (iter)->GetAtomInMolecule() * (iter)->GetAtomicNum() *
                (iter)->GetAverageNucleonWeight();
    }

    sumCharge_ = aux1;
    sumNucleons_ = aux2;
    ZA_ = aux1 / aux2;
    molDensity_ = massDensity_ * NA / aux2;
    MM_ = aux3 / aux2;

    // TODO: Compare to Bremsstrahlung::CalculateScatteringX0; just one (this or
    // the Brems-thing) is needed Calculation of the radiation length
    aux1 = 0;
    aux2 = 0;

    for (std::vector<Component>::iterator iter =
             components_.begin();
         iter != components_.end(); ++iter) {
        aux1 += X0_inv((iter)->GetNucCharge(), (iter)->GetAtomicNum()) *
                ((iter)->GetAtomInMolecule() * (iter)->GetAtomicNum());
        aux2 += (iter)->GetAtomInMolecule() * (iter)->GetAtomicNum();
    }

    radiationLength_ = aux2 / aux1;
    //radiationLength_ /= massDensity_; //New radiation length is in grammage
}

// ------------------------------------------------------------------------- //
// Getter
// ------------------------------------------------------------------------- //


size_t Medium::GetHash() const noexcept
{
    size_t hash_digest = 0;
    hash_combine(hash_digest, numComponents_, sumCharge_, ZA_, I_, C_, a_, m_,
        X0_, X1_, d0_, massDensity_, molDensity_, radiationLength_, MM_,
        sumNucleons_);
    for (auto component : components_)
        hash_combine(hash_digest, component.GetHash());
    return hash_digest;
}

Medium Medium::GetMediumForHash(size_t hash) {
    auto it = Medium::medium_map->find(hash);
    if (it != Medium::medium_map->end())
        return it->second;
    throw std::invalid_argument("Medium for given hash not registered");
}

/******************************************************************************
 *                              Different Media                               *
 ******************************************************************************/

PDG2001::Water::Water()
    : Medium(
        "water",
        75.0,     // I
        -3.5017,  // C
        0.09116,  // a
        3.4773,   // m
        0.2400,   // X0
        2.8004,   // X1
        0,        // d0
        1.000,    // massDensitiy
        {Components::Hydrogen(2),
         Components::Oxygen()}) {}

PDG2020::Water::Water()
    : Medium(
          "water",
          79.7,     // I
          -3.5017,  // C
          0.09116,  // a
          3.4773,   // m
          0.2400,   // X0
          2.8004,   // X1
          0,        // d0
          1.000,    // massDensitiy
          {Components::Hydrogen(2),
           Components::Oxygen()}) {}

PDG2001::Ice::Ice()
    : Medium(
          "ice",
          75.0,     // I
          -3.5017,  // C
          0.09116,  // a
          3.4773,   // m
          0.2400,   // X0
          2.8004,   // X1
          0,        // d0
          0.917,    // massDensitiy
          {Components::Hydrogen(2),
           Components::Oxygen()}) {}

PDG2020::Ice::Ice()
    : Medium(
          "ice",
          79.7,     // I
          -3.5873,  // C
          0.09116,  // a
          3.4773,   // m
          0.2586,   // X0
          2.8190,   // X1
          0,        // d0
          0.918,    // massDensitiy
          {Components::Hydrogen(2),
           Components::Oxygen()}) {}

Salt::Salt()
    : Medium("salt",
             // Calculated by ESTAR detabase
             // (it could be 185 eV by the method of reference below)
             175.3,    // I
             -4.5041,  // C
             0.1632,   // a
             3,        // m
             0.2,      // X0
             3.0,      // X1
             0,        // d0
             2.323,    // Solid halite density
             {Components::Sodium(), Components::Chlorine()}) {}

CalciumCarbonate::CalciumCarbonate()
    : Medium(
          "calciumcarbonate",
          136.4,    // I
          -3.7738,  // C
          0.08301,  // a
          3.4120,   // m
          0.0492,   // X0
          3.0549,   // X1
          0,        // d0
          2.650,    // massDensity
          {Components::Calcium(),
           Components::Carbon(),
           Components::Oxygen(3)}) {
}

StandardRock::StandardRock()
    : Medium("standardrock",
             136.4,    // I
             -3.7738,  // C
             0.08301,  // a
             3.4120,   // m
             0.0492,   // X0
             3.0549,   // X1
             0,        // d0
             2.650,    // massDensity
             {Components::StandardRock()}) {}

FrejusRock::FrejusRock()
    : Medium("frejusrock",
             149.0,   // I
             -5.053,  // C
             0.078,   // a
             3.645,   // m
             0.288,   // X0
             3.196,   // X1
             0,       // d0
             2.740,   // massDensity
             {Components::FrejusRock()}) {}

Iron::Iron()
    : Medium("iron",
             286.0,    // I
             -4.2911,  // C
             0.14680,  // a
             2.9632,   // m
             -0.0012,  // X0
             3.1531,   // X1
             0.12,     // d0
             7.874,    // massDensity
             {Components::Iron()}) {
}

Hydrogen::Hydrogen()
    : Medium("hydrogen",
             21.8,     // I
             -3.0977,  // C
             0.13483,  // a
             5.6249,   // m
             0.4400,   // X0
             1.8856,   // X1
             0,        // d0
             0.07080,  // massDensity
             {Components::Hydrogen()}) {}

Lead::Lead()
    : Medium("lead",
             823.0,    // I
             -6.2018,  // C
             0.09359,  // a
             3.1608,   // m
             0.3776,   // X0
             3.8073,   // X1
             0.14,     // d0
             11.350,   // massDensity
             {Components::Lead()}) {
}

Copper::Copper()
    : Medium(
          "copper",
          322.0,    // I
          -4.4190,  // C
          0.14339,  // a
          2.9044,   // m
          -0.0254,  // X0
          3.2792,   // X1
          0.08,     // d0
          8.960,    // massDensity
          {Components::Copper()}) {}

Uranium::Uranium()
    : Medium(
          "uranium",
          890.0,    // I
          -5.8694,  // C
          0.19677,  // a
          2.8171,   // m
          0.2260,   // X0
          3.3721,   // X1
          0.14,     // d0
          18.950,   // massDensity
          {Components::Uranium()}) {
}

Air::Air()
    : Medium("air",
             85.7,      // I
             -10.5961,  // C
             0.10914,   // a
             3.3994,    // m
             1.7418,    // X0
             4.2759,    // X1
             0,         // d0
             1.205e-3,  // dry, 1 atm massDensity
             {Components::Nitrogen(2 * 78.1),
              Components::Oxygen(2 * 21.0),
              Components::Argon(0.9)}) {}

Paraffin::Paraffin()
    : Medium(
          "paraffin",
          55.9,     // I
          -2.9551,  // C
          0.1209,   // a
          3.4288,   // m
          0.1289,   // X0
          2.5084,   // X1
          0,        // d0
          0.93,     // massDensity
          {Components::Carbon(25.0),
           Components::Oxygen(52.0)}) {}

AntaresWater::AntaresWater()
    : Medium("antareswater",
             75.0,     // I
             -3.5017,  // C
             0.09116,  // a
             3.4773,   // m
             0.2400,   // X0
             2.8004,   // X1
             0,        // d0
             1.03975,  // massDensity
             {Components::Hydrogen(2.0),
              Components::Oxygen(1.00884),
              Components::Sodium(0.00943),
              Components::Potassium(0.000209),
              Components::Magnesium(0.001087),
              Components::Calcium(0.000209),
              Components::Chlorine(0.01106),
              Components::Sulfur(0.00582)}) {
    /*
     * initialize ANTARES water
     * Sea water (Mediterranean Sea, ANTARES place)
     * ==========================================================================
     * WATER DENSITY CHANGES WITH THE DEPTH FROM 1.0291 g/cm^3 AT SURFACE
     * UP TO 1.0404 g/cm^3 AT THE SEA BED
     * (ANTARES-Site/2000-001 and references therein)
     *
     * The error which is caused by this simplified approach (average value for
     * density) does not exceed 0.5% (much less, in fact) that is comparable
     *with an error which comes from uncertainties with the muon cross-sections.
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
CascadiaBasinWater::CascadiaBasinWater()
    : Medium(
        "cascadiabasinwater",
        75.,       // I
        -3.5017,   // C
        0.09116,   // a
        3.4773,    // m
        0.2400,    // X0
        2.8004,    // X1
        0.,        // d0
        1.0400322, // massDensity
        {
            Components::Hydrogen(2.0),
            Components::Oxygen(1.0021),
            Components::Sodium(8.7174e-3),
            Components::Magnesium(9.8180e-4),
            Components::Calcium(1.9113e-4),
            Components::Potassium(1.8975e-4),
            Components::Chlorine(1.0147e-2),
            Components::Sulfur(5.2487e-4)
        })
{}

LiquidArgon::LiquidArgon()
    : Medium("liquidargon",
             188.0,     // I
             -5.2146,   // C
             0.19559,   // a
             3.0000,    // m
             0.2000,    // X0
             3.0000,    // X1
             0,         // d0
             1.396,     // massDensity
             {Components::Argon()}) {}

/******************************************************************************
 *                        private Helper Funcitons                             *
 ******************************************************************************/

double Medium::X0_inv(unsigned int Z, double M) {
    double a_sq = 0.;
    double fZ = 0.;
    double Lrad = 0.;
    double Lrad_dash = 0.;

    a_sq = ALPHA * ALPHA * Z * Z;
    fZ = a_sq * (1. / (1. + a_sq) + 0.20206 - 0.0369 * a_sq +
                 0.0083 * a_sq * a_sq - 0.002 * a_sq * a_sq * a_sq);

    // Get radiation logarithm from Tsai (Rev.Mod.Phys.46)
    if (Z > 4)  // Elements Z>4
    {
        // Thomas-Fermi-Moliere model (Tsai eq. B21)
        Lrad = std::log(184.15 * std::pow(Z, -1. / 3.));
        Lrad_dash = std::log(1194. * std::pow(Z, -2. / 3.));
    }
    // Tsai table B2
    else if (Z == 1)  // Hydrogen
    {
        Lrad = 5.31;
        Lrad_dash = 6.144;
    } else if (Z == 2)  // Helium
    {
        Lrad = 4.79;
        Lrad_dash = 5.621;
    } else if (Z == 3)  // Lithium
    {
        Lrad = 4.74;
        Lrad_dash = 5.805;
    } else if (Z == 4)  // Beryllium
    {
        Lrad = 4.71;
        Lrad_dash = 5.924;
    }
    // Bremsstrahlung complete screening case
    return 4. * ALPHA * RE * RE * NA / M *
           (Z * Z * (Lrad - fZ) + Z * Lrad_dash);
}

namespace PROPOSAL {
const std::string GetMediumName( int enumVal ) { return Medium_Name[enumVal]; }
} // namespace PROPOSAL
