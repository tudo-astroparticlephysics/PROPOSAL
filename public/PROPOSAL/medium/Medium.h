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


namespace PROPOSAL {

namespace Components
{
    class Component;
}


/******************************************************************************
*                                   Medium                                    *
******************************************************************************/

class Medium
{
    public:
    /*!
     * initialize medium by its name andby the
     * multiplicative density correction factor
     * \param   w       medium to create
     * \param   rho     multiplicative density correction factor
     */
    // TODO(mario): Doc string Thu 2017/08/03
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
           double massDensity);
    Medium(const Medium&);

    ///@brief Crush this Medium.
    virtual ~Medium();

    void swap(Medium& medium);
    virtual Medium* clone() const = 0; // Prototyping/Virtual constructor idiom (used for deep copies)

    // Operators
    Medium& operator=(const Medium&);
    bool operator==(const Medium& medium) const;
    bool operator!=(const Medium& medium) const;
    friend std::ostream& operator<<(std::ostream& os, Medium const& medium);

    void init(); // Needed when not using c++11

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

    double X0_inv(unsigned int Z, double M);
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
class MediumCopyable: virtual public Base
{
    public:
    virtual Base* clone() const
    {
        return new Derived(static_cast<Derived> (*this));
    }

    static  Base* create(double density_correction = 1.0)
    {
        return new Derived(density_correction);
    }

    Derived& operator=(const Medium& medium)
    {
        if (this != &medium)
        {
            const Derived* med = static_cast<Derived const&>(&medium);  // Throws on bad cast.

            Derived tmp(*med);
            swap(tmp);
        }

        return *this;
    }
};

// ----------------------------------------------------------------------------
/// @brief Implement Medium Water
// ----------------------------------------------------------------------------
class Water : public MediumCopyable<Medium, Water>
{
    public:
    Water(double rho = 1.0);
    Water(const Medium& medium)
        : Medium(medium)
    {
    }
    virtual ~Water() {}
};

// ----------------------------------------------------------------------------
/// @brief Implement Medium Ice
// ----------------------------------------------------------------------------
class Ice : public MediumCopyable<Medium, Ice>
{
    public:
    Ice(double rho = 1.0);
    Ice(const Medium& medium)
        : Medium(medium)
    {
    }
    virtual ~Ice() {}
};

// ----------------------------------------------------------------------------
/// @brief Implement Medium Salt
// ----------------------------------------------------------------------------
class Salt : public MediumCopyable<Medium, Salt>
{
    public:
    Salt(double rho = 1.0);
    Salt(const Medium& medium)
        : Medium(medium)
    {
    }
    virtual ~Salt() {}
};

// ----------------------------------------------------------------------------
/// @brief Implement Medium CalciumCarbonate (CaCO3)
// ----------------------------------------------------------------------------
class CalciumCarbonate : public MediumCopyable<Medium, CalciumCarbonate>
{
    public:
    CalciumCarbonate(double rho = 1.0);
    CalciumCarbonate(const Medium& medium)
        : Medium(medium)
    {
    }
    virtual ~CalciumCarbonate() {}
};

// ----------------------------------------------------------------------------
/// @brief Implement Medium StandardRock
// ----------------------------------------------------------------------------
class StandardRock : public MediumCopyable<Medium, StandardRock>
{
    public:
    StandardRock(double rho = 1.0);
    StandardRock(const Medium& medium)
        : Medium(medium)
    {
    }
    virtual ~StandardRock() {}
};

// ----------------------------------------------------------------------------
/// @brief Implement Medium FrejusRock
// ----------------------------------------------------------------------------
class FrejusRock : public MediumCopyable<Medium, FrejusRock>
{
    public:
    FrejusRock(double rho = 1.0);
    FrejusRock(const Medium& medium)
        : Medium(medium)
    {
    }
    virtual ~FrejusRock() {}
};

// ----------------------------------------------------------------------------
/// @brief Implement Medium Iron
// ----------------------------------------------------------------------------
class Iron : public MediumCopyable<Medium, Iron>
{
    public:
    Iron(double rho = 1.0);
    Iron(const Medium& medium)
        : Medium(medium)
    {
    }
    virtual ~Iron() {}
};

// ----------------------------------------------------------------------------
/// @brief Implement Medium Hydrogen
// ----------------------------------------------------------------------------
class Hydrogen : public MediumCopyable<Medium, Hydrogen>
{
    public:
    Hydrogen(double rho = 1.0);
    Hydrogen(const Medium& medium)
        : Medium(medium)
    {
    }
    virtual ~Hydrogen() {}
};

// ----------------------------------------------------------------------------
/// @brief Implement Medium Lead
// ----------------------------------------------------------------------------
class Lead : public MediumCopyable<Medium, Lead>
{
    public:
    Lead(double rho = 1.0);
    Lead(const Medium& medium)
        : Medium(medium)
    {
    }
    virtual ~Lead() {}
};

// ----------------------------------------------------------------------------
/// @brief Implement Medium Copper
// ----------------------------------------------------------------------------
class Copper : public MediumCopyable<Medium, Copper>
{
    public:
    Copper(double rho = 1.0);
    Copper(const Medium& medium)
        : Medium(medium)
    {
    }
    virtual ~Copper() {}
};

// ----------------------------------------------------------------------------
/// @brief Implement Medium Uranium
// ----------------------------------------------------------------------------
class Uranium : public MediumCopyable<Medium, Uranium>
{
    public:
    Uranium(double rho = 1.0);
    Uranium(const Medium& medium)
        : Medium(medium)
    {
    }
    virtual ~Uranium() {}
};

// ----------------------------------------------------------------------------
/// @brief Implement Medium Air
// ----------------------------------------------------------------------------
class Air : public MediumCopyable<Medium, Air>
{
    public:
    Air(double rho = 1.0);
    Air(const Medium& medium)
        : Medium(medium)
    {
    }
    virtual ~Air() {}
};

// ----------------------------------------------------------------------------
/// @brief Implement Medium Paraffin
// ----------------------------------------------------------------------------
class Paraffin : public MediumCopyable<Medium, Paraffin>
{
    public:
    Paraffin(double rho = 1.0);
    Paraffin(const Medium& medium)
        : Medium(medium)
    {
    }
    virtual ~Paraffin() {}
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
class AntaresWater : public MediumCopyable<Medium, AntaresWater>
{
    public:
    AntaresWater(double rho = 1.0);
    AntaresWater(const Medium& medium)
        : Medium(medium)
    {
    }
    virtual ~AntaresWater() {}
};

} // namespace PROPOSAL
