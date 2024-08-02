
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/

#pragma once

#include "PROPOSAL/medium/Components.h"
#include <array>
#include <cmath>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#define MEDIUM_DEF(cls)                                                        \
    class cls : public Medium {                                                \
    public:                                                                    \
        cls();                                                                 \
    };

namespace PROPOSAL {

/******************************************************************************
 *                                   Medium                                    *
 ******************************************************************************/

struct IonizationParameters {
    const double I; // mean excitation energy [eV]
    const double C;
    const double a; //< ionization formula constants
    const double m;
    const double X0;
    const double X1; //< ionization formula constants (continued)
    const double d0; // sternheimer et al density effect parameters
public:
    // IonizationParameters() = default;

    IonizationParameters(
        double I, double C, double a, double m, double X0, double X1, double d0)
        : I(I)
        , C(C)
        , a(a)
        , m(m)
        , X0(X0)
        , X1(X1)
        , d0(d0)
    {
    }

    static IonizationParameters approx(
        const std::vector<Component>& comp, double density)
    {
        auto I = calculate_mean_ionization_energy(comp);
        auto ZA = calculate_mean_ZA(comp);
        auto plasma_energy = calculate_plasma_energy(ZA, density);
        auto cBar = calculate_cBar(plasma_energy, I);
        auto m = 3.f;
        auto x1 = calculate_x1(I);
        auto x0 = calculate_x0(I, cBar);
        auto a = calculate_a(cBar, x0, x1);
        auto d0 = 0.f;
        return IonizationParameters(I, -cBar, a, m, x0, x1, d0);
    }

    static std::vector<double> calculate_weights(
        const std::vector<Component>& comp)
    {
        auto weights = std::vector<double>(comp.size());
        for (const auto& c : comp)
            weights.push_back(c.GetAtomInMolecule() * c.GetNucCharge());
        auto sum = std::accumulate(weights.begin(), weights.end(), 0.f);
        for (auto& w : weights)
            w /= sum;
        return weights;
    }

    static double calculate_mean_ionization_energy(
        const std::vector<Component>& comp)
    {
        auto aux1 = 0.f;
        auto aux2 = 0.f;
        for (const auto& c : comp) {
            aux1 += c.GetAtomInMolecule() * c.GetAtomicNum()
                * std::log(c.GetIonizationEnergy());
            aux2 += c.GetAtomInMolecule() * c.GetAtomicNum();
        }

        return std::exp(aux1 / aux2);
    }

    static double calculate_mean_ZA(const std::vector<Component>& comp)
    {
        auto A = std::accumulate(comp.cbegin(), comp.cend(), 0.f,
            [](double sum, const Component& c) {
                return sum + c.GetAtomInMolecule() * c.GetAtomicNum();
            });
        auto Z = std::accumulate(comp.cbegin(), comp.cend(), 0.f,
            [](double sum, const Component& c) {
                return sum + c.GetAtomInMolecule() * c.GetNucCharge();
            });
        return Z / A;
    }

    static double calculate_plasma_energy(double ZA, double density)
    {
        return 28.816 * std::sqrt(ZA * density);
    }

    static double calculate_cBar(double plasma_energy, double I)
    {
        return 2 * std::log(I / plasma_energy) + 1;
    }

    static double calculate_x1(double I)
    {
        if (I < 100)
            return 2.0;
        else
            return 3.0;
    }

    static double calculate_x0(double I, double cBar)
    {
        if (I < 100) {
            if (cBar < 3.681)
                return 0.2;
            return 0.326 * cBar - 1.0;
        }
        if (cBar < 5.215)
            return 0.2;
        return 0.326 * cBar - 1.5;
    }

    static double calculate_a(double cBar, double x0, double x1)
    {
        return (cBar - 2 * std::log(10) * x0) / std::pow(x1 - x0, 3);
    }
};

class Medium {
    IonizationParameters ioniz;

public:
    // class Builder;

public:
    // Medium() = default;
    Medium(std::string name, double I, double C, double a, double m, double X0,
        double X1, double d0, double massDensity, std::vector<Component>);

    Medium(std::string name, double massDensity, std::vector<Component> comp);

    // Operators
    bool operator==(const Medium& medium) const;
    bool operator!=(const Medium& medium) const;
    friend std::ostream& operator<<(std::ostream& os, Medium const& medium);

    // ----------------------------------------------------------------- //
    // Getter & Setter
    // ----------------------------------------------------------------- //

    // Getter
    int GetNumComponents() const { return components_.size(); }
    std::vector<Component> GetComponents() const { return components_; }
    double GetSumCharge() const { return sumCharge_; }
    double GetZA() const { return sumCharge_ / sumNucleons_; }
    double GetI() const { return ioniz.I; }
    double GetC() const { return ioniz.C; }
    double GetA() const { return ioniz.a; }
    double GetM() const { return ioniz.m; }
    double GetX0() const { return ioniz.X0; }
    double GetX1() const { return ioniz.X1; }
    double GetD0() const { return ioniz.d0; }
    double GetMassDensity() const { return massDensity_; }
    double GetRadiationLength() const { return radiationLength_; }
    double GetMolDensity() const;
    std::string GetName() const { return name_; }
    double GetMM() const { return MM_; }
    double GetSumNucleons() const { return sumNucleons_; }
    size_t GetHash() const noexcept;

    static Medium GetMediumForHash(size_t);

protected:
    // Methods
    void init();
    // double X0_inv(unsigned int Z, double M);

    // Protected member
    std::string name_;

    int numComponents_; ///< number of components

    double sumCharge_;   ///< sum of charges of all nuclei
    double sumNucleons_; ///< sum of nucleons of all nuclei

    double ZA_; // <Z/A>
    // double I_;           // mean excitation energy [eV]
    // double C_, a_;       //< ionization formula constants
    // double m_, X0_, X1_; //< ionization formula constants (continued)
    // double d0_;          // sternheimer et al density effect parameters

    double massDensity_; ///< mass density (20° C, 1 atm) [g/cm3]
    // double molDensity_;      ///< molecule density [number/cm3]
    double radiationLength_; ///< radiation length [cm]

    double MM_; ///< average all-component nucleon weight

    std::vector<Component> components_; ///< Components of Medium

private:
    static std::unique_ptr<std::map<size_t, Medium>> medium_map;
};

namespace PDG2001 {
    MEDIUM_DEF(Water)
    MEDIUM_DEF(Ice)
}
namespace PDG2020 {
    MEDIUM_DEF(Water)
    MEDIUM_DEF(Ice)
}
using namespace PDG2020;

MEDIUM_DEF(Salt)
MEDIUM_DEF(CalciumCarbonate)
MEDIUM_DEF(StandardRock)
MEDIUM_DEF(FrejusRock)
MEDIUM_DEF(Iron)
MEDIUM_DEF(Hydrogen)
MEDIUM_DEF(Lead)
MEDIUM_DEF(Copper)
MEDIUM_DEF(Uranium)
MEDIUM_DEF(Paraffin)
MEDIUM_DEF(Air)
MEDIUM_DEF(LiquidArgon)

// Plagioclase 39%
MEDIUM_DEF(Albit)
MEDIUM_DEF(Anorthit)

// Alkali feldspars 12%
MEDIUM_DEF(Sanidine)
MEDIUM_DEF(Anorthoclase)

// Quartz 12%
MEDIUM_DEF(Quartz)

// Pyroxene 11%
MEDIUM_DEF(Enstatite)
MEDIUM_DEF(Ferrosilit)

// Amphibole 5%
MEDIUM_DEF(Anthophyllit)

// Mica / Glimmer 5 %
MEDIUM_DEF(Muscovite)
MEDIUM_DEF(Phlogopite)

// Olivine Group 3%
MEDIUM_DEF(Liebenbergite)

// clay mineral 4.5%
MEDIUM_DEF(Kaolinite)

// carbonate annd nitrate 2%
MEDIUM_DEF(Calcit)
MEDIUM_DEF(Dolomite)

// Magnetit 1.5%
MEDIUM_DEF(Magnetit)

// sulphate mineral
MEDIUM_DEF(Gypsum)

MEDIUM_DEF(Cobaltite)
MEDIUM_DEF(Linneit)

// #<{(|
// * initialize ANTARES water
// * Sea water (Mediterranean Sea, ANTARES place)
// * ==========================================================================
// * WATER DENSITY CHANGES WITH THE DEPTH FROM 1.0291 g/cm^3 AT SURFACE
// * UP TO 1.0404 g/cm^3 AT THE SEA BED
// * (ANTARES-Site/2000-001 and references therein)
// *
// * The error which is caused by this simplified approach (average value
// for
// * density) does not exceed 0.5% (much less, in fact) that is comparable
// with
// *  an error which comes from uncertainties with the muon cross-sections.
// *==========================================================================
// |)}>#
MEDIUM_DEF(AntaresWater)

/// @brief Cascadia Basin seawater
///
/// Cascadia Basin is located in the Northeast Pacific Ocean at a depth of
/// 2.66 km close to the coast of Victoria, Canada. Infrastructure for
/// deep-sea experiments is provided by Ocean Networks Canada (ONC).
/// Therefore, it is an interesting site for future neutrino telescopes.
///
/// The chemical composition of the Cacadia Basin seawater is based on the
/// six most important (by reference salinity) materials dissolved in
/// Standard Seawater, according to F. J. Millero et al., Deep Sea Research
/// Part I: Oceanographic Research Papers 55.1 (2008), pp. 50-72. These
/// materials are sodium, magnesium, calcium, potassium, chlorine, and
/// sulfate ions.
///
/// The project https://kkrings.github.io/seawater/ is used to calculate the
/// chemical composition of the Cascadia Basin seawater. The practical
/// salinity of 34.6288 psu and the mass density of 1.0400322 g/cm3 are
/// obtained from ONC measurements.
MEDIUM_DEF(CascadiaBasinWater)

} // namespace PROPOSAL
#undef MEDIUM_DEF

namespace PROPOSAL {
enum Medium_Type {
    WATER,
    WATERPDG2001,
    WATERPDG2020,
    ICE,
    ICEPDG2001,
    ICEPDG2020,
    SALT,
    STANDARDROCK,
    FREJUSROCK,
    IRON,
    HYDROGEN,
    LEAD,
    COPPER,
    URANIUM,
    AIR,
    PARAFFIN,
    ANTARESWATER,
    CASCADIABASINWATER,
    LIQUIDARGON,
    ALBIT,
    ANORTHIT,
    SANIDINE,
    ANORTHOCLASE,
    QUARTZ,
    ENSTATITE,
    FERROSILIT,
    ANTHOPHYLLIT,
    MUSCOVITE,
    PHLOGOPITE,
    LIEBENBERGITE,
    KAOLINITE,
    CALCIT,
    DOLOMITE,
    MAGNETIT,
    GYPSUM,
    COBALTITE,
    LINNEIT,
};
} // namespace PROPOSAL

namespace PROPOSAL {
static const std::array<std::string, 19> Medium_Name = {
    "water",
    "waterpdg2001",
    "waterpdg2020",
    "ice",
    "icepdg2001",
    "icepdg2020",
    "salt",
    "standardrock",
    "frejusrock",
    "iron",
    "hydrogen",
    "lead",
    "copper",
    "uranium",
    "air",
    "paraffin",
    "antareswater",
    "cascadiabasinwater",
    "liquidargon",
};
} // namespace PROPOSAL

/* namespace PROPOSAL { */
/* const std::string GetMediumName( int enumVal ); */
/* } // namespace PROPOSAL */
