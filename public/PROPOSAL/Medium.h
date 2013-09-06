/*! \file   Medium.h
*   \brief  Header file for definition of the medium class object.
*
*   For more details see the class documentation.
*
*   \date   21.06.2010
*   \author Jan-Hendrik Koehne
*/

#ifndef Medium_H
#define Medium_H

#include <vector>
#include <string>
#include <boost/bind.hpp>
#include <boost/function.hpp>

class Medium
{

protected:

    int numCompontents_;                	///< number of components
    std::vector<double> nucCharge_;         ///< nucleus charge
    std::vector<double> atomicNum_;         ///< atomic number
    std::vector<double> atomInMolecule_;    ///< number of atoms in a molecule
    double sumCharge_;                  	///< sum of charges of all nuclei

    double ZA_;                         	///< <Z/A>
    double I_;                          	///< ionization potential [eV]
    double C1_, C_, a_;                 	///< ionization formula constants
    double m_, X0_, X1_, d0_;           	///< ionization formula constants (continued)
    double r_;                          	///< refraction index

    std::vector<double> logConstant_;       ///< radiation logarithm constant B
    std::vector<double> bPrime_;            ///< radiation logarithm constant bPrime

    double rho_;                        	///< multiplicative density correction factor
    double massDensity_;                	///< mass density [g/cm3]
    double molDensity_;                 	///< molecule density [number/cm3]
    double radiationLength_;                ///< radiation length [cm]

    std::vector<double> M_;                 ///< average nucleon weight in a nucleus [MeV]
    std::vector<std::string> elementName_;  ///< element name
    std::string name_;                  	///< medium name

    double ecut_;                       	///< cutoff energy [MeV]
    double vcut_;                       	///< relative cutoff energy
    double vCut_;                       	///< relative cutoff energy - call setCut(E) to set this

    std::vector<double> mN_;                ///< Woods-Saxon potential factor
    double MM_;                         	///< average all-component nucleon weight
    double sumNucleons_;                	///< sum of nucleons of all nuclei

    double r0_;

//----------------------------------------------------------------------------//

    // Memberfunctions

    /*!
     * set the value of radiation logarithm constant B
     *
     * \param   i   nucleon number
     * \return  value of radiation logarithm constant B
     */

    void SetLogConstant(int i);
//----------------------------------------------------------------------------//
    /*!
     * set the value of radiation logarithm constant bPrime
     *
     * \param   i   nucleon number
     * \return  value of radiation logarithm constant bPrime
     */

    void SetBPrime(int i);

//----------------------------------------------------------------------------//
    /*!
     * initialization of arrays
     *
     * \param   i       number of components
     */

    void Inita(int i);

//----------------------------------------------------------------------------//
    /*!
     *initialization code, common to all media
     */

    void Initr();

//----------------------------------------------------------------------------//
    /*!
     * initialize water
     *
     */

    void InitWater();
//----------------------------------------------------------------------------//
    /*!
     * initialize ice
     */

    void InitIce();
//----------------------------------------------------------------------------//
    /*!
     * initialize salt (added by Ped)
     */

    void InitSalt();
//----------------------------------------------------------------------------//
    /*!
     * initialize standard rock
     */

    void InitStandardrock();
//----------------------------------------------------------------------------//
    /*!
     * initialize Frejus rock
     */

    void InitFrejusrock();
//----------------------------------------------------------------------------//
    /*!
     * initialize iron
     */

    void InitIron();
//----------------------------------------------------------------------------//
    /*!
     * initialize hydrogen
     */

    void InitHydrogen();
//----------------------------------------------------------------------------//
    /*!
     * initialize lead
     */

    void InitLead();
//----------------------------------------------------------------------------//
    /*!
     * initialize copper
     */

    void InitCopper();
//----------------------------------------------------------------------------//
    /*!
     * initialize uranium
     */

    void InitUranium();
//----------------------------------------------------------------------------//
    /*!
     * initialize air
     */

    void InitAir();
//----------------------------------------------------------------------------//
    /*!
     * initialize mineral oil or paraffin CH3(CH2)~23CH3 (added by Ped)
     */

    void InitParaffin();
//----------------------------------------------------------------------------//
    /*!
     * initialize ANTARES water
     * Sea water (Mediterranean Sea, ANTARES place) \n
     * ==========================================================================\n
     * WATER DENSITY CHANGES WITH THE DEPTH FROM 1.0291 g/cm^3 AT SURFACE
     * UP TO 1.0404 g/cm^3 AT THE SEA BED
     * (ANTARES-Site/2000-001 and references therein) \n

     * The error which is caused by this simplified approach (average value for
     * density) does not exceed 0.5% (much less, in fact) that is comparable with
     * an error which comes from uncertainties with the muon cross-sections.
     *  ==========================================================================\n
     * \n
     * added by Claudine Colnard, \n
     * Institute Nikhef, The Netherlands, \n
     * ANTARES collaboration.
     */

    void InitAntaresWater();
//----------------------------------------------------------------------------//


    /*!
     * Woods-Saxon potential calculation - function to integrate
     *
     * \param   r
     * \return  value of the Woods-Saxon potential
     */

    double FunctionToIntegral(double r);

//----------------------------------------------------------------------------//

public:

    // constructors

    Medium();
    Medium(const Medium&);
    Medium& operator=(const Medium&);
    bool operator==(const Medium &medium) const;
    bool operator!=(const Medium &medium) const;
    friend std::ostream& operator<<(std::ostream& os, Medium const& medium);


//----------------------------------------------------------------------------//

    /*!
     * initialize medium by its name andby the
     * multiplicative density correction factor
     * \param   w       medium to create
     * \param   rho     multiplicative density correction factor
     */


    Medium(std::string w, double rho);

//----------------------------------------------------------------------------//

   void swap(Medium &medium);

//----------------------------------------------------------------------------//
    // Getter

    int GetNumComponents() const
    {
        return numCompontents_;
    }
//----------------------------------------------------------------------------//

    std::vector<double> GetNucCharge() const
    {
        return nucCharge_;
    }

 //----------------------------------------------------------------------------//

    std::vector<double> GetAtomicNum() const
    {
        return atomicNum_;
    }
//----------------------------------------------------------------------------//

    std::vector<double> GetAtomInMolecule() const
    {
        return atomInMolecule_;
    }
//----------------------------------------------------------------------------//

    double GetSumCharge() const
    {
        return sumCharge_;
    }
//----------------------------------------------------------------------------//
    double GetZA() const
    {
        return ZA_;
    }
//----------------------------------------------------------------------------//
    double GetI() const
    {
        return I_;
    }
//----------------------------------------------------------------------------//
    double GetC1() const
    {
        return C1_;
    }
//----------------------------------------------------------------------------//
    double GetC() const
    {
        return C_;
    }
//----------------------------------------------------------------------------//
    double GetA() const
    {
        return a_;
    }
//----------------------------------------------------------------------------//
    double GetM() const
    {
        return m_;
    }
//----------------------------------------------------------------------------//
    double GetX0() const
    {
        return X0_;
    }
//----------------------------------------------------------------------------//
    double GetX1() const
    {
        return X1_;
    }
//----------------------------------------------------------------------------//
    double GetD0() const
    {
        return d0_;
    }
//----------------------------------------------------------------------------//
    double GetR() const
    {
        return r_;
    }
//----------------------------------------------------------------------------//
    std::vector<double> GetLogConstant() const
    {
        return logConstant_;
    }
//----------------------------------------------------------------------------//
    std::vector<double> GetBPrime() const
    {
        return bPrime_;
    }
//----------------------------------------------------------------------------//
    double GetRho() const
    {
        return rho_;
    }
//----------------------------------------------------------------------------//
    double GetMassDensity() const
    {
        return massDensity_;
    }
//----------------------------------------------------------------------------//
    double GetRadiationLength() const
    {
        return radiationLength_;
    }
//----------------------------------------------------------------------------//
    double GetMolDensity() const
    {
        return molDensity_;
    }
//----------------------------------------------------------------------------//
    std::vector<double> GetAverageNucleonWeight() const
    {
        return M_;
    }
//----------------------------------------------------------------------------//
    std::vector<std::string> GetElementName() const
    {
        return elementName_;
    }
//----------------------------------------------------------------------------//
    std::string GetElementName(int i) const
    {
        return elementName_.at(i);
    }
//----------------------------------------------------------------------------//
    std::string GetName() const
    {
        return name_;
    }
//----------------------------------------------------------------------------//
    std::vector<double> GetMN() const
    {
        return mN_;
    }
//----------------------------------------------------------------------------//
    double GetMM() const
    {
        return MM_;
    }
//----------------------------------------------------------------------------//
    double GetSumNucleons() const
    {
        return sumNucleons_;
    }
//----------------------------------------------------------------------------//
    double GetR0() const
    {
        return r0_;
    }
//----------------------------------------------------------------------------//

    // Setter
    void SetNumComponents(int numCompontents);
    void SetNucCharge(std::vector<double> nucCharge);
    void SetAtomicNum(std::vector<double> atomicNum);
    void SetAtomInMolecule(std::vector<double> atomInMolecule);
    void SetSumCharge(double sumCharge);
    void SetZA(double ZA);
    void SetI(double I);
    void SetC1(double C1);
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
    void SetElementName(std::vector<std::string> E);
    void SetName(std::string name);
    void SetMN(std::vector<double> mN);
    void SetMM(double MM);
    void SetSumNucleons(double sumNucleons);
    void SetR0(double r0);

//----------------------------------------------------------------------------//
    // destructors

    ///@brief Crush this Medium.
    ~Medium();

};



#endif //Medium_H
