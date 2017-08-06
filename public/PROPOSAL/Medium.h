/*! \file   Medium.h
*   \brief  Header file for definition of the medium class object.
*
*   For more details see the class documentation.
*
*   \date   21.06.2010
*   \author Jan-Hendrik Koehne
*/

#pragma once

#include <map>
#include <string>
#include <vector>

#include "PROPOSAL/Components.h"

namespace PROPOSAL {

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

    void swap(Medium& medium);
    virtual Medium* clone() const = 0; // Prototyping/Virtual constructor idiom (used for deep copies)

    // Operators
    Medium& operator=(const Medium&);
    bool operator==(const Medium& medium) const;
    bool operator!=(const Medium& medium) const;
    friend std::ostream& operator<<(std::ostream& os, Medium const& medium);

    void init(); // Needed when not using c++11

    ///@brief Crush this Medium.
    virtual ~Medium();

    // ----------------------------------------------------------------- //
    // Getter & Setter
    // ----------------------------------------------------------------- //

    // Getter
    int GetNumComponents() const { return numComponents_; }

    std::vector<Components::Component*>& GetComponents() { return components_; }

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

    double GetRho() const { return rho_; }

    double GetMassDensity() const { return massDensity_; }

    double GetRadiationLength() const { return radiationLength_; }

    double GetMolDensity() const { return molDensity_; }

    std::string GetName() const { return name_; }

    double GetMM() const { return MM_; }

    double GetSumNucleons() const { return sumNucleons_; }

    double GetR0() const { return r0_; }

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
    void SetR0(double r0);

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

    double r0_;
};

/******************************************************************************
*                               MediumCopyable                                *
******************************************************************************/


// ----------------------------------------------------------------------------
/// @brief Template class for Medium
///
/// Provides a template to the assignment operator
// ----------------------------------------------------------------------------
template<typename T>
class MediumCopyable: public Medium
{
    public:
    MediumCopyable(std::string name,
                   double rho,
                   double I,
                   double C,
                   double a,
                   double m,
                   double X0,
                   double X1,
                   double d0,
                   double massDensity)
        : Medium(name, rho, I, C, a, m, X0, X1, d0, massDensity)
    {
    }

    T& operator=(const Medium& medium)
    {
        if (this != &medium)
        {
            const T* med = dynamic_cast<T const&>(&medium);  // Throws on bad cast.

            T tmp(*med);
            swap(tmp);
        }

        return *this;
    }
};

// #<{(|
// * initialize water
// |)}>#
// class Water : public Medium
class Water : public MediumCopyable<Water>
{
    public:
    Water(double rho = 1.0);
    Water(const Water& medium)
        : MediumCopyable(medium)
    {
    }
    virtual ~Water() {}

    Medium* clone() const { return new Water(*this); };
    static Medium* create() { return new Water(); };
};

// #<{(|
// * initialize ice
// |)}>#
class Ice : public MediumCopyable<Ice>
{
    public:
    Ice(double rho = 1.0);
    Ice(const Ice& medium)
        : MediumCopyable(medium)
    {
    }
    virtual ~Ice() {}

    Medium* clone() const { return new Ice(*this); };
    static Medium* create() { return new Ice(); };
};

// #<{(|
// * initialize salt (added by Ped)
// |)}>#
class Salt : public MediumCopyable<Salt>
{
    public:
    Salt(double rho = 1.0);
    Salt(const Salt& medium)
        : MediumCopyable(medium)
    {
    }
    virtual ~Salt() {}

    Medium* clone() const { return new Salt(*this); };
    static Medium* create() { return new Salt(); };
};

// #<{(|
// * initialize standard rock
// |)}>#
class StandardRock : public MediumCopyable<StandardRock>
{
    public:
    StandardRock(double rho = 1.0);
    StandardRock(const StandardRock& medium)
        : MediumCopyable(medium)
    {
    }
    virtual ~StandardRock() {}

    Medium* clone() const { return new StandardRock(*this); };
    static Medium* create() { return new StandardRock(); };
};

// #<{(|
// * initialize Frejus rock
// |)}>#
class FrejusRock : public MediumCopyable<FrejusRock>
{
    public:
    FrejusRock(double rho = 1.0);
    FrejusRock(const FrejusRock& medium)
        : MediumCopyable(medium)
    {
    }
    virtual ~FrejusRock() {}

    Medium* clone() const { return new FrejusRock(*this); };
    static Medium* create() { return new FrejusRock(); };
};

// #<{(|
// * initialize iron
// |)}>#
class Iron : public MediumCopyable<Iron>
{
    public:
    Iron(double rho = 1.0);
    Iron(const Iron& medium)
        : MediumCopyable(medium)
    {
    }
    virtual ~Iron() {}

    Medium* clone() const { return new Iron(*this); };
    static Medium* create() { return new Iron(); };
};

// #<{(|
// * initialize hydrogen
// |)}>#
class Hydrogen : public MediumCopyable<Hydrogen>
{
    public:
    Hydrogen(double rho = 1.0);
    Hydrogen(const Hydrogen& medium)
        : MediumCopyable(medium)
    {
    }
    virtual ~Hydrogen() {}

    Medium* clone() const { return new Hydrogen(*this); };
    static Medium* create() { return new Hydrogen(); };
};

// #<{(|
// * initialize lead
// |)}>#
class Lead : public MediumCopyable<Lead>
{
    public:
    Lead(double rho = 1.0);
    Lead(const Lead& medium)
        : MediumCopyable(medium)
    {
    }
    virtual ~Lead() {}

    Medium* clone() const { return new Lead(*this); };
    static Medium* create() { return new Lead(); };
};

// #<{(|
// * initialize copper
// |)}>#
class Copper : public MediumCopyable<Copper>
{
    public:
    Copper(double rho = 1.0);
    Copper(const Copper& medium)
        : MediumCopyable(medium)
    {
    }
    virtual ~Copper() {}

    Medium* clone() const { return new Copper(*this); };
    static Medium* create() { return new Copper(); };
};

// #<{(|
// * initialize uranium
// |)}>#
class Uranium : public MediumCopyable<Uranium>
{
    public:
    Uranium(double rho = 1.0);
    Uranium(const Uranium& medium)
        : MediumCopyable(medium)
    {
    }
    virtual ~Uranium() {}

    Medium* clone() const { return new Uranium(*this); };
    static Medium* create() { return new Uranium(); };
};

// #<{(|
// * initialize air
// |)}>#
class Air : public MediumCopyable<Air>
{
    public:
    Air(double rho = 1.0);
    Air(const Air& medium)
        : MediumCopyable(medium)
    {
    }
    virtual ~Air() {}

    Medium* clone() const { return new Air(*this); };
    static Medium* create() { return new Air(); };
};

// #<{(|
// * initialize mineral oil or paraffin CH3(CH2)~23CH3 (added by Ped)
// |)}>#
class Paraffin : public MediumCopyable<Paraffin>
{
    public:
    Paraffin(double rho = 1.0);
    Paraffin(const Paraffin& medium)
        : MediumCopyable(medium)
    {
    }
    virtual ~Paraffin() {}

    Medium* clone() const { return new Paraffin(*this); };
    static Medium* create() { return new Paraffin(); };
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
class AntaresWater : public MediumCopyable<AntaresWater>
{
    public:
    AntaresWater(double rho = 1.0);
    AntaresWater(const AntaresWater& medium)
        : MediumCopyable(medium)
    {
    }
    virtual ~AntaresWater() {}

    Medium* clone() const { return new AntaresWater(*this); };
    static Medium* create() { return new AntaresWater(); };
};

/******************************************************************************
*                               Medium Factory                                *
******************************************************************************/

class MediumFactory
{
    public:
    virtual ~MediumFactory() { medium_map.clear(); }

    void Register(const std::string& name, Medium* (*)(void));
    Medium* CreateMedium(const std::string&);

    static MediumFactory* Get()
    {
        static MediumFactory instance;
        return &instance;
    }

    private:
    MediumFactory();
    std::map<std::string, Medium* (*)(void)> medium_map;
};

} // namespace PROPOSAL
