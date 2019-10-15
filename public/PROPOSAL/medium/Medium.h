
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

#include <string>
#include <vector>
#include <memory>

#include "PROPOSAL/medium/Components.h"

#define MEDIUM_DEF(cls)                                                                                                \
    class cls : public Medium                                                                                          \
    {                                                                                                                  \
    public:                                                                                                            \
        cls(double rho = 1.0);                                                                                         \
        cls(const Medium& medium)                                                                                      \
            : Medium(medium)                                                                                           \
        {                                                                                                              \
        }                                                                                                              \
        virtual ~cls() {}                                                                                              \
                                                                                                                       \
        virtual Medium* clone() const { return new cls(*this); }                                                       \
        static Medium* create(double density_correction = 1.0) { return new cls(density_correction); }                 \
    };

namespace PROPOSAL {

/******************************************************************************
 *                                   Medium                                    *
 ******************************************************************************/

class Medium
{
public:
    // class Builder;

public:
    Medium() {}
    Medium(std::string name,
           double rho,
           double I,
           double C,
           double a,
           double m,
           double X0,
           double X1,
           double d0,
           double massDensity,
           std::vector<std::shared_ptr<Components::Component>>);
    Medium(const Medium&);
    virtual Medium* clone() const { return new Medium(*this); };

    ///@brief Crush this Medium.
    virtual ~Medium();

    void swap(Medium& medium);

    // Operators
    Medium& operator=(const Medium&);
    bool operator==(const Medium& medium) const;
    bool operator!=(const Medium& medium) const;
    friend std::ostream& operator<<(std::ostream& os, Medium const& medium);

    // ----------------------------------------------------------------- //
    // Getter & Setter
    // ----------------------------------------------------------------- //

    // Getter
    int GetNumComponents() const { return numComponents_; }
    const std::vector<Components::Component*>& GetComponents() const { return components_; }
    double GetSumCharge() const { return sumCharge_; }
    double GetZA() const { return ZA_; }
    double GetI() const { return I_; }
    double GetC() const { return C_; }
    double GetA() const { return a_; }
    double GetM() const { return m_; }
    double GetX0() const { return X0_; }
    double GetX1() const { return X1_; }
    double GetD0() const { return d0_; }
    double GetR() const { return r_; }
    double GetDensityCorrection() const { return rho_; }
    double GetMassDensity() const { return massDensity_; }
    double GetRadiationLength() const { return radiationLength_; }
    double GetMolDensity() const { return molDensity_; }
    std::string GetName() const { return name_; }
    double GetMM() const { return MM_; }
    double GetSumNucleons() const { return sumNucleons_; }

    // Setter
    void SetNumComponents(int numComponents);
    void SetSumCharge(double sumCharge);
    void SetZA(double ZA);
    void SetI(double I);
    void SetC(double C);
    void SetA(double a);
    void SetM(double m);
    void SetX0(double X0);
    void SetX1(double X1);
    void SetD0(double d0);
    void SetR(double r);
    void SetRho(double rho);
    void SetMassDensity(double massDensity);
    void SetMolDensity(double molDensity);
    void SetAverageNucleonWeight(std::vector<double> M);
    void SetComponents(std::vector<std::shared_ptr<Components::Component>> components);
    void SetMM(double MM);
    void SetSumNucleons(double sumNucleons);

protected:
    // Methods
    void init();
    double X0_inv(unsigned int Z, double M);

    // Protected member
    std::string name_;

    int numComponents_;                              ///< number of components
    std::vector<Components::Component*> components_; ///< Components of Medium

    double sumCharge_; ///< sum of charges of all nuclei

    double ZA_;               ///< <Z/A>
    double I_;                ///< ionization potential [eV]
    double C_, a_;            ///< ionization formula constants
    double m_, X0_, X1_, d0_; ///< ionization formula constants (continued)
    double r_;                ///< refraction index

    double rho_;             ///< multiplicative density correction factor
    double massDensity_;     ///< mass density [g/cm3]
    double molDensity_;      ///< molecule density [number/cm3]
    double radiationLength_; ///< radiation length [cm]

    double ecut_; ///< cutoff energy [MeV]
    double vcut_; ///< relative cutoff energy
    double vCut_; ///< relative cutoff energy - call setCut(E) to set this

    double MM_;          ///< average all-component nucleon weight
    double sumNucleons_; ///< sum of nucleons of all nuclei
};

MEDIUM_DEF(Water)
MEDIUM_DEF(Ice)
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

class Air : public Medium
{
public:
    static const double fraction_N;
    static const double fraction_O;
    static const double fraction_Ar;
    static const double fraction_sum;

public:
    Air(double rho = 1.0);
    Air(const Medium& medium)
        : Medium(medium)
    {
    }
    virtual ~Air() {}

    virtual Medium* clone() const { return new Air(*this); }
    static Medium* create(double density_correction = 1.0) { return new Air(density_correction); }
};

// #<{(|
// * initialize ANTARES water
// * Sea water (Mediterranean Sea, ANTARES place)
// * ==========================================================================
// * WATER DENSITY CHANGES WITH THE DEPTH FROM 1.0291 g/cm^3 AT SURFACE
// * UP TO 1.0404 g/cm^3 AT THE SEA BED
// * (ANTARES-Site/2000-001 and references therein)
// *
// * The error which is caused by this simplified approach (average value for
// * density) does not exceed 0.5% (much less, in fact) that is comparable with
// *  an error which comes from uncertainties with the muon cross-sections.
// *==========================================================================
// |)}>#
MEDIUM_DEF(AntaresWater)

/// @brief Cascadia Basin seawater
///
/// Cascadia Basin is located in the Northeast Pacific Ocean at a depth of
/// 2.66 km close to the coast of Victoria, Canada. Infrastructure for deep-sea
/// experiments is provided by Ocean Networks Canada (ONC). Therefore, it is an
/// interesting site for future neutrino telescopes.
///
/// The chemical composition of the Cacadia Basin seawater is based on the six
/// most important (by reference salinity) materials dissolved in Standard
/// Seawater, according to F. J. Millero et al., Deep Sea Research Part I:
/// Oceanographic Research Papers 55.1 (2008), pp. 50-72. These materials are
/// sodium, magnesium, calcium, potassium, chlorine, and sulfate ions.
///
/// The project https://kkrings.github.io/seawater/ is used to calculate the
/// chemical composition of the Cascadia Basin seawater. The practical salinity
/// of 34.6288 psu and the mass density of 1.0400322 g/cm3 are obtained from
/// ONC measurements.
MEDIUM_DEF(CascadiaBasinWater)

} // namespace PROPOSAL

#undef MEDIUM_DEF
