/*! \file   Medium.h
*   \brief  Header file for definition of the medium class object.
*
*   For more details see the class documentation.
*
*   \date   21.06.2010
*   \author Jan-Hendrik Koehne
*/
#pragma once

#ifndef Medium_H
#define Medium_H

#include <vector>
#include <string>

namespace PROPOSAL
{

namespace Components
{

class Component
{
    public:
    Component(std::string name,
              double charge,
              double atomicNum,
              double atomInMolecule);
    virtual ~Component() = 0;

    // Getter
    std::string GetName() const { return name_; }
    double GetNucCharge() const { return nucCharge_; }
    double GetAtomicNum() const { return atomicNum_; }
    double GetAtomInMolecule() const { return atomInMolecule_; }
    double GetLogConstat() const { return logConstant_; }
    double GetBPrime() const { return bPrime_; }
    double GetAverageNucleonWeight() const { return M_; }
    double GetWoodSaxon() const { return mN_; }
    double GetR0() const { return r0_; }

    protected:
    /*!
     * set the value of radiation logarithm constant B
     *
     * \param   i   nucleon number
     * \return  value of radiation logarithm constant B
     */
    void SetLogConstant();

    /*!
     * set the value of radiation logarithm constant bPrime
     *
     * \param   i   nucleon number
     * \return  value of radiation logarithm constant bPrime
     */
    void SetBPrime();

    /*!
     * Woods-Saxon potential calculation - function to integrate
     *
     * \param   r
     * \return  value of the Woods-Saxon potential
     */
    double FunctionToIntegral(double r);

    // Passed to constructor
    std::string name_;
    double nucCharge_;      ///< nucleus charge
    double atomicNum_;      ///< atomic number
    double atomInMolecule_; ///< number of atoms in a molecule

    // Calculated in constructor
    double logConstant_; ///< radiation logarithm constant B
    double bPrime_;      ///< radiation logarithm constant bPrime
    double M_;           ///< average nucleon weight in a nucleus [MeV]
    double mN_;          ///< Woods-Saxon potential factor
    double r0_;          // //TODO(mario): Must realy be stored? Thu 2017/08/03
};

class Oxygen: public Component
{
    public:
        Oxygen(double atomInMolecule = 1.0);
        virtual ~Oxygen();
};

class Hydrogen: public Component
{
    public:
        Hydrogen(double atomInMolecule = 2.0);
        virtual ~Hydrogen();
};

class Natrium: public Component
{
    public:
        Natrium(double atomInMolecule = 1.0);
        virtual ~Natrium();
};

class Chloride: public Component
{
    public:
        Chloride(double atomInMolecule = 1.0);
        virtual ~Chloride();
};

class Standard_Rock: public Component
{
    public:
        Standard_Rock(double atomInMolecule = 1.0);
        virtual ~Standard_Rock();
};

class Frejus_Rock: public Component
{
    public:
        Frejus_Rock(double atomInMolecule = 1.0);
        virtual ~Frejus_Rock();
};

class Iron: public Component
{
    public:
        Iron(double atomInMolecule = 1.0);
        virtual ~Iron();
};

class Lead: public Component
{
    public:
        Lead(double atomInMolecule = 1.0);
        virtual ~Lead();
};

class Copper: public Component
{
    public:
        Copper(double atomInMolecule = 1.0);
        virtual ~Copper();
};

class Uranium: public Component
{
    public:
        Uranium(double atomInMolecule = 1.0);
        virtual ~Uranium();
};

class Nitrogen: public Component
{
    public:
        Nitrogen(double atomInMolecule = 1.0);
        virtual ~Nitrogen();
};

class Arsenic: public Component
{
    public:
        Arsenic(double atomInMolecule = 1.0);
        virtual ~Arsenic();
};

class Carbon: public Component
{
    public:
        Carbon(double atomInMolecule = 1.0);
        virtual ~Carbon();
};

class Potassium: public Component
{
    public:
        Potassium(double atomInMolecule = 1.0);
        virtual ~Potassium();
};

class Magnesium: public Component
{
    public:
        Magnesium(double atomInMolecule = 1.0);
        virtual ~Magnesium();
};

class Calcium: public Component
{
    public:
        Calcium(double atomInMolecule = 1.0);
        virtual ~Calcium();
};

class Sulfur: public Component
{
    public:
        Sulfur(double atomInMolecule = 1.0);
        virtual ~Sulfur();
};

}

/******************************************************************************
*                                   Medium                                    *
******************************************************************************/


class Medium
{
    public:

        // Medium();
        Medium(const Medium&);
        /*!
         * initialize medium by its name andby the
         * multiplicative density correction factor
         * \param   w       medium to create
         * \param   rho     multiplicative density correction factor
         */
        //TODO(mario): Doc string Thu 2017/08/03
        Medium(std::string name,
               double rho,
               std::vector<Components::Component*>,
               double I,
               double C,
               double a,
               double m,
               double X0,
               double X1,
               double d0,
               double massDensity);

        void swap(Medium &medium);

        // Operators
        // Medium& operator=(const Medium&);
        bool operator==(const Medium &medium) const;
        bool operator!=(const Medium &medium) const;
        friend std::ostream& operator<<(std::ostream& os, Medium const& medium);

        ///@brief Crush this Medium.
        virtual ~Medium() = 0;  // Pure virtual to ensure no pure Medium will be created

        // ----------------------------------------------------------------- //
        // Getter & Setter
        // ----------------------------------------------------------------- //


        // Getter
        int GetNumComponents() const
        {
            return numComponents_;
        }

        double GetSumCharge() const
        {
            return sumCharge_;
        }

        double GetZA() const
        {
            return ZA_;
        }

        double GetI() const
        {
            return I_;
        }

        double GetC() const
        {
            return C_;
        }

        double GetA() const
        {
            return a_;
        }

        double GetM() const
        {
            return m_;
        }

        double GetX0() const
        {
            return X0_;
        }

        double GetX1() const
        {
            return X1_;
        }

        double GetD0() const
        {
            return d0_;
        }

        double GetR() const
        {
            return r_;
        }

        double GetRho() const
        {
            return rho_;
        }

        double GetMassDensity() const
        {
            return massDensity_;
        }

        double GetRadiationLength() const
        {
            return radiationLength_;
        }

        double GetMolDensity() const
        {
            return molDensity_;
        }

        std::string GetName() const
        {
            return name_;
        }

        double GetMM() const
        {
            return MM_;
        }

        double GetSumNucleons() const
        {
            return sumNucleons_;
        }

        double GetR0() const
        {
            return r0_;
        }


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
        // void SetName(std::string name);
        void SetComponents(std::vector<Components::Component*>);
        void SetMM(double MM);
        void SetSumNucleons(double sumNucleons);
        void SetR0(double r0);

    protected:
        // ----------------------------------------------------------------- //
        // Protected Memeber
        // ----------------------------------------------------------------- //

        std::string name_;

        int numComponents_;                            ///< number of components
        std::vector<Components::Component*> components_; ///< Components of Medium
        // std::vector<double> nucCharge_;         ///< nucleus charge
        // std::vector<double> atomicNum_;         ///< atomic number
        // std::vector<double> atomInMolecule_;    ///< number of atoms in a
        // molecule
        double sumCharge_; ///< sum of charges of all nuclei

        double ZA_;                         	///< <Z/A>
        double I_;                          	///< ionization potential [eV]
        double C_, a_;                          ///< ionization formula constants
        double m_, X0_, X1_, d0_;           	///< ionization formula constants (continued)
        double r_;                          	///< refraction index

        double rho_;                        	///< multiplicative density correction factor
        double massDensity_;                	///< mass density [g/cm3]
        double molDensity_;                 	///< molecule density [number/cm3]
        double radiationLength_;                ///< radiation length [cm]

        double ecut_;                       	///< cutoff energy [MeV]
        double vcut_;                       	///< relative cutoff energy
        double vCut_;                       	///< relative cutoff energy - call setCut(E) to set this

        double MM_;                         	///< average all-component nucleon weight
        double sumNucleons_;                	///< sum of nucleons of all nuclei

        double r0_;



    // //----------------------------------------------------------------------------//
    //     #<{(|!
    //      * initialize water
    //      *
    //      |)}>#
    //
    //     void InitWater();
    // //----------------------------------------------------------------------------//
    //     #<{(|!
    //      * initialize ice
    //      |)}>#
    //
    //     void InitIce();
    // //----------------------------------------------------------------------------//
    //     #<{(|!
    //      * initialize salt (added by Ped)
    //      |)}>#
    //
    //     void InitSalt();
    // //----------------------------------------------------------------------------//
    //     #<{(|!
    //      * initialize standard rock
    //      |)}>#
    //
    //     void InitStandardrock();
    // //----------------------------------------------------------------------------//
    //     #<{(|!
    //      * initialize Frejus rock
    //      |)}>#
    //
    //     void InitFrejusrock();
    // //----------------------------------------------------------------------------//
    //     #<{(|!
    //      * initialize iron
    //      |)}>#
    //
    //     void InitIron();
    // //----------------------------------------------------------------------------//
    //     #<{(|!
    //      * initialize hydrogen
    //      |)}>#
    //
    //     void InitHydrogen();
    // //----------------------------------------------------------------------------//
    //     #<{(|!
    //      * initialize lead
    //      |)}>#
    //
    //     void InitLead();
    // //----------------------------------------------------------------------------//
    //     #<{(|!
    //      * initialize copper
    //      |)}>#
    //
    //     void InitCopper();
    // //----------------------------------------------------------------------------//
    //     #<{(|!
    //      * initialize uranium
    //      |)}>#
    //
    //     void InitUranium();
    // //----------------------------------------------------------------------------//
    //     #<{(|!
    //      * initialize air
    //      |)}>#
    //
    //     void InitAir();
    // //----------------------------------------------------------------------------//
    //     #<{(|!
    //      * initialize mineral oil or paraffin CH3(CH2)~23CH3 (added by Ped)
    //      |)}>#
    //
    //     void InitParaffin();
    // //----------------------------------------------------------------------------//
    //     #<{(|!
    //      * initialize ANTARES water
    //      * Sea water (Mediterranean Sea, ANTARES place) \n
    //      * ==========================================================================\n
    //      * WATER DENSITY CHANGES WITH THE DEPTH FROM 1.0291 g/cm^3 AT SURFACE
    //      * UP TO 1.0404 g/cm^3 AT THE SEA BED
    //      * (ANTARES-Site/2000-001 and references therein) \n
    //
    //      * The error which is caused by this simplified approach (average value for
    //      * density) does not exceed 0.5% (much less, in fact) that is comparable with
    //      * an error which comes from uncertainties with the muon cross-sections.
    //      *  ==========================================================================\n
    //      * \n
    //      * added by Claudine Colnard, \n
    //      * Institute Nikhef, The Netherlands, \n
    //      * ANTARES collaboration.
    //      |)}>#
    //
    //     void InitAntaresWater();
    //----------------------------------------------------------------------------//
};


class Water: public Medium
{
    public:
        Water(double rho);
        virtual ~Water();
};

}


#endif //Medium_H
