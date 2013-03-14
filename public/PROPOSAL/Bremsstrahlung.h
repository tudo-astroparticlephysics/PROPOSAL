#include <vector>
#include "CrossSections.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/Particle.h"

class Bremsstrahlung: public CrossSections
{
private:

    bool        init_;
    double      eLpm_;
    double      xo_;

    bool        lorenz_;
    double      lorenzCut_;  	/// in [MeV] // - set to 1.e6 in Constructor


    std::vector<Integral*>   integrals_;

//----------------------------------------------------------------------------//
    //Memberfunctions

    // The following functions define the different
    // parametrization for the Bremsstrahlung Cross Section
    /*!
    \param v relativ energy loss
    \param i the nucleon in the medium on which the bremsstahlung occur
    \return  Calculates \f$ a_{1} \f$ see function Sel in this class
    */
    double KelnerKakoulinPetrukhinParametrization(double v, int i);

//----------------------------------------------------------------------------//
    /*!
    \param v relativ energy loss
    \param i the nucleon in the medium on which the bremsstahlung occur
    \return Calculates \f$ a_{1} \f$ see function Sel in this class
    */

    double AndreevBezrukovBugaevParametrization(double v, int i);

//----------------------------------------------------------------------------//
    /*!
    \param v relativ energy loss
    \param i the nucleon in the medium on which the bremsstahlung occur
    \return Calculates \f$ a_{1} \f$ see function Sel in this class
    */
    double PetrukhinShestakovParametrization(double v, int i);

    //----------------------------------------------------------------------------//
    /*!
    \param v relativ energy loss
    \param i the nucleon in the medium on which the bremsstahlung occur
    \return Calculates \f$ a_{1} \f$ see function Sel in this class
    */
    double CompleteScreeningCase(double v, int i);

public:
//----------------------------------------------------------------------------//

    Bremsstrahlung();
    Bremsstrahlung(const Bremsstrahlung&);
    Bremsstrahlung& operator=(const Bremsstrahlung&);

//----------------------------------------------------------------------------//

    void SetIntegralLimits();

//----------------------------------------------------------------------------//

    double CalculatedEdx();

//----------------------------------------------------------------------------//

    double CalculatedNdx();


//----------------------------------------------------------------------------//

    double CalculatedNdx(double rnd);

//----------------------------------------------------------------------------//

    double CalculateStochasticLoss();

//----------------------------------------------------------------------------//

    void EnableStochasticInerpolation();

//----------------------------------------------------------------------------//

    void EnableContinuousInerpolation();

//----------------------------------------------------------------------------//

    double FunctionToContinuousIntegral(double variable);

//----------------------------------------------------------------------------//

    double FunctionToStochasticalIntegral(double variable);

//----------------------------------------------------------------------------//

    /*!
    this is what the Elastic Bremsstrahlung Cross Section (EBCS) is equal to
    units are [1/cm] since the multiplication by No*n is done here.
    Corrections for excitations of the nucleus and deep inelastic
    excitations of separate nucleons are
    included (positive term dn/Z), as well as the contribution of
    the mu-diagrams to the inelastic
    bremsstrahlung on the electrons (non-zero only for allowed
    energies of photon after electron recoil).
    four different parametrizations are enumerated, the final result is:
    \f[ \sigma_{el}= \rho_{mol}N(z^2)^2\Big[a_1\Big] \f] and \f$ a_{1} \f$
    depends on the chosen parametrization
    \param i the nucleon in the medium on which the bremsstahlung occur
    \param v relativ energy loss
    \return Elastic Bremsstrahlung Cross Section [1/cm]
    */

    double ElasticBremsstrahlungCrossSection(double v, int i);

//----------------------------------------------------------------------------//
    ~Bremsstrahlung(){}

};
