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
#include "PROPOSAL/Integral.h"

#include "PhysicsModel.h"

class Medium: public PhysicsModel
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
    std::vector<double> M_;                 ///< average nucleon weight in a nucleus [MeV]
    std::vector<std::string> E_;            ///< element name
    std::string name_;                  	///< medium name

    double ecut_;                       	///< cutoff energy [MeV]
    double vcut_;                       	///< relative cutoff energy
    double vCut_;                       	///< relative cutoff energy - call setCut(E) to set this

    Integral integral_;                 	// Needed for calculation of the
    std::vector<double> mN_;                ///< Woods-Saxon potential factor
    double MM_;                         	///< average all-component nucleon weight
    double sumNucleons_;                	///< sum of nucleons of all nuclei

    double r0_;

//----------------------------------------------------------------------------//

public:

    // constructors

    /// @brief  Create an Medium.
    Medium();
//----------------------------------------------------------------------------//

    /*!
     * initialize medium by its name and proposed values of energy cut (ecut)
     * and fractional energy cut (vcut).
     * The cuts ecut and vcut to be used must satisfy \f$e_{cut}\geq 0\f$ and
     * \f$0 \leq v_{cut} \leq 1\f$; If both values satisfy these
     * inequalities, the lower from \f$ E \cdot v_{cut} \f$ and \f$ e_{cut}\f$
     * will be used. If only one value satisfies the inequalities,
     * only that one value is used. If both values are outside these intervals,
     * \f$v_{cut}=1\f$ is assumed.
     *
     * \param   w       medium to create
     * \param   ecut    stochastic calculation when energy > ecut
     * \param   vcut    stochastic calculation when e-loss > energy*vcut
     */


    Medium(std::string w, double ecut, double vcut, double rho);

//----------------------------------------------------------------------------//

    // Memberfunctions

    /*!
     * initialization of arrays
     *
     * \param   i       number of components
     */

    void inita(int i);

//----------------------------------------------------------------------------//
    /*!
     *initialization code, common to all media
     */

    void initr();

//----------------------------------------------------------------------------//
    /*!
     * set the value of radiation logarithm constant B
     *
     * \param   i   nucleon number
     * \return  value of radiation logarithm constant B
     */

    double elB(int i);
//----------------------------------------------------------------------------//
    /*!
     * set the value of radiation logarithm constant bPrime
     *
     * \param   i   nucleon number
     * \return  value of radiation logarithm constant bPrime
     */

    double elP(int i);
//----------------------------------------------------------------------------//
    /*!
     * initialize water
     *
     */

    void initWater();
//----------------------------------------------------------------------------//
    /*!
     * initialize ice
     */

    void initIce();
//----------------------------------------------------------------------------//
    /*!
     * initialize salt (added by Ped)
     */

    void initSalt();
//----------------------------------------------------------------------------//
    /*!
     * initialize standard rock
     */

    void initStandardrock();
//----------------------------------------------------------------------------//
    /*!
     * initialize Frejus rock
     */

    void initFrejusrock();
//----------------------------------------------------------------------------//
    /*!
     * initialize iron
     */

    void initIron();
//----------------------------------------------------------------------------//
    /*!
     * initialize hydrogen
     */

    void initHydrogen();
//----------------------------------------------------------------------------//
    /*!
     * initialize lead
     */

    void initLead();
//----------------------------------------------------------------------------//
    /*!
     * initialize uranium
     */

    void initUranium();
//----------------------------------------------------------------------------//
    /*!
     * initialize air
     */

    void initAir();
//----------------------------------------------------------------------------//
    /*!
     * initialize mineral oil or paraffin CH3(CH2)~23CH3 (added by Ped)
     */

    void initParaffin();
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

    void initAntaresWater();
//----------------------------------------------------------------------------//

    /*!
     * return the value of the fractional energy cut,
     * as described in the help for the constructor
     *
     * \param   E   Energy
     * \return  value of the fractional energy cut
     */

    double vCut(double E);
//----------------------------------------------------------------------------//
    /*!
     * Woods-Saxon potential calculation - interface to Integral
     *
     * \param   r
     * \return  value of the Woods-Saxon potential
     */

    double function(double r);
//----------------------------------------------------------------------------//
    // Getter

    int get_numCompontents() const
    {
        return numCompontents_;
    }

    std::vector<double> get_NucCharge() const
    {
        return nucCharge_;
    }

    std::vector<double> get_atomicNum() const
    {
        return atomicNum_;
    }

    std::vector<double> get_atomInMolecule() const
    {
        return atomInMolecule_;
    }

    double get_sumCharge() const
    {
        return sumCharge_;
    }

    double get_ZA() const
    {
        return ZA_;
    }

    double get_I() const
    {
        return I_;
    }

    double get_C1() const
    {
        return C1_;
    }

    double get_C() const
    {
        return C_;
    }

    double get_a() const
    {
        return a_;
    }

    double get_m() const
    {
        return m_;
    }

    double get_X0() const
    {
        return X0_;
    }

    double get_X1() const
    {
        return X1_;
    }

    double get_d0() const
    {
        return d0_;
    }

    double get_r() const
    {
        return r_;
    }

    std::vector<double> get_logConstant() const
    {
        return logConstant_;
    }

    std::vector<double> get_bPrime() const
    {
        return bPrime_;
    }

    double get_rho() const
    {
        return rho_;
    }

    double get_massDensity() const
    {
        return massDensity_;
    }

    double get_molDensity() const
    {
        return molDensity_;
    }

    std::vector<double> get_M() const
    {
        return M_;
    }

    std::vector<std::string> get_E() const
    {
        return E_;
    }

    std::string get_E(int i) const
    {
        return E_[i];
    }

    std::string get_name() const
    {
        return name_;
    }

    double get_ecut() const
    {
        return ecut_;
    }

    double get_vcut() const
    {
        return vcut_;
    }

    double get_vCut() const
    {
        return vCut_;
    }

    Integral get_integral() const
    {
        return integral_;
    }

    std::vector<double> get_mN() const
    {
        return mN_;
    }

    double get_MM() const
    {
        return MM_;
    }

    double get_sumNucleons() const
    {
        return sumNucleons_;
    }

    double get_r0() const
    {
        return r0_;
    }

//----------------------------------------------------------------------------//
    // Setter <--- declared but not defined yet

    void set_numCompontents(int numCompontents);
    void set_NucCharge(std::vector<double> NucCharge);
    void set_atomicNum(std::vector<double> atomicNum);
    void set_atomInMolecule(std::vector<double> atomInMolecule);
    void set_sumCharge(double sumCharge);
    void set_ZA(double ZA);
    void set_I(double I);
    void set_C1(double C1);
    void set_C(double C);
    void set_a(double a);
    void set_m(double m);
    void set_X0(double X0);
    void set_X1(double X1);
    void set_d0(double d0);
    void set_r(double r);
    void set_logConstant(std::vector<double> logConstant);
    void set_bPrime(std::vector<double> bPrime);
    void set_rho(double rho);
    void set_massDensity(double massDensity);
    void set_molDensity(double molDensity);
    void set_M(std::vector<double> M);
    void set_E(std::vector<std::string> E);
    void set_name(std::string name);
    void set_ecut(double ecut);
    void set_vcut(double vcut);
    void set_vCut(double vCut);
    void set_integral(Integral integral);
    void set_mN(std::vector<double> mN);
    void set_MM(double MM);
    void set_sumNucleons(double sumNucleons);
    void set_r0(double r0);

//----------------------------------------------------------------------------//
    // destructors

    ///@brief Crush this Medium.
    ~Medium();


};



#endif //Medium_H
