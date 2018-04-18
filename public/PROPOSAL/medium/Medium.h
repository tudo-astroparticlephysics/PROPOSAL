/*! \file   Medium.h
 *   \brief  Header file for definition of the medium class object.
 *
 *   For more details see the class documentation.
 *
 *   \date   21.06.2010
 *   \author Jan-Hendrik Koehne
 */

#pragma once

#include <string>
#include <vector>

#include "PROPOSAL/medium/Components.h"

#define MEDIUM_DEF(cls)                                                                                                \
    class cls : public MediumCopyable<Medium, cls>                                                                     \
    {                                                                                                                  \
    public:                                                                                                            \
        cls(double rho = 1.0);                                                                                         \
        cls(const Medium& medium)                                                                                      \
            : Medium(medium)                                                                                           \
        {                                                                                                              \
        }                                                                                                              \
        virtual ~cls() {}                                                                                              \
    };

namespace PROPOSAL {

/******************************************************************************
 *                                   Medium                                    *
 ******************************************************************************/

class Medium
{
public:
    class Builder;

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
           const std::vector<Components::Component*>&);
    Medium(const Medium&);
    virtual Medium* clone() const
    {
        return new Medium(*this);
    }; // Prototyping/Virtual constructor idiom (used for deep copies)

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
    void SetComponents(std::vector<Components::Component*>);
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

/******************************************************************************
 *                                  Builder                                    *
 ******************************************************************************/

class Medium::Builder
{
public:
    Builder();

    // --------------------------------------------------------------------- //
    // Setter
    // --------------------------------------------------------------------- //

    Builder& SetName(const std::string& var)
    {
        this->name_ = var;
        return *this;
    }
    Builder& SetRho(double var)
    {
        rho_ = var;
        return *this;
    }
    Builder& SetI(double var)
    {
        I_ = var;
        return *this;
    }
    Builder& SetC(double var)
    {
        C_ = var;
        return *this;
    }
    Builder& SetA(double var)
    {
        a_ = var;
        return *this;
    }
    Builder& SetM(double var)
    {
        m_ = var;
        return *this;
    }
    Builder& SetX0(double var)
    {
        X0_ = var;
        return *this;
    }
    Builder& SetX1(double var)
    {
        X1_ = var;
        return *this;
    }
    Builder& SetD0(double var)
    {
        d0_ = var;
        return *this;
    }
    Builder& SetMassDensity(double var)
    {
        massDensity_ = var;
        return *this;
    }
    Builder& addComponent(const Components::Component& var)
    {
        components_.push_back(var.clone());
        return *this;
    }

    Builder& SetMedium(const Medium& var)
    {
        name_        = var.GetName();
        rho_         = var.GetDensityCorrection();
        I_           = var.GetI();
        C_           = var.GetC();
        a_           = var.GetA();
        m_           = var.GetM();
        X0_          = var.GetX0();
        X1_          = var.GetX1();
        d0_          = var.GetD0();
        massDensity_ = var.GetMassDensity();

        components_ = var.GetComponents();
        return *this;
    }

    Medium build() { return Medium(name_, rho_, I_, C_, a_, m_, X0_, X1_, d0_, massDensity_, components_); }

private:
    std::string name_;
    double rho_;
    double I_;
    double C_;
    double a_;
    double m_;
    double X0_;
    double X1_;
    double d0_;
    double massDensity_;

    std::vector<Components::Component*> components_;
};

/******************************************************************************
 *                               MediumCopyable                                *
 ******************************************************************************/

// ----------------------------------------------------------------------------
/// @brief Template class for Medium
///
/// Provides a template to the assignment operator, clone & create
// ----------------------------------------------------------------------------
template<class Base, class Derived>
class MediumCopyable : virtual public Base
{
public:
    virtual Base* clone() const { return new Derived(static_cast<Derived>(*this)); }

    static Base* create(double density_correction = 1.0) { return new Derived(density_correction); }

    Derived& operator=(const Medium& medium)
    {
        if (this != &medium)
        {
            const Derived* med = static_cast<Derived const&>(&medium); // Throws on bad cast.

            Derived tmp(*med);
            swap(tmp);
        }

        return *this;
    }
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

class Air : public MediumCopyable<Medium, Air>
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

} // namespace PROPOSAL

#undef MEDIUM_DEF
