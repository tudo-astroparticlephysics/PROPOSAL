/*
 * Epairproduction.h
 *
 *  Created on: 29.06.2010
 *      Author: koehne
 */

#ifndef EPAIRPRODUCTION_H_
#define EPAIRPRODUCTION_H_


#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Interpolate.h"

#include "CrossSections.h"

/*! \class Epairproduction Epairproduction.h "Epairproduction.h"
    \brief class contains functions necessary for calculation of e+e- pair production losses
 */


class EpairContinuous;
class EpairStochastic;

class Epairproduction: public CrossSections
{


protected:



    double vMax, vUp, vMin;
    Epairproduction *epair;
    Integral *integral_;
    int i;
    double v;
    bool init;
    double eLpm;
    bool jt_;

//----------------------------------------------------------------------------//



public:



    Interpolate* interpolateJ_;
    EpairContinuous *continuous_;
    EpairStochastic *stochastic_;

//----------------------------------------------------------------------------//

    //Constructors

    /*!
    creates internal references to p and m, to be called from subclasses
    */

    Epairproduction(Epairproduction *cros);

//----------------------------------------------------------------------------//

    /*!
    initializes subclasses and creates internal references to p and m
    */

    Epairproduction(CrossSections *cros);

//----------------------------------------------------------------------------//
    /*!
    call before using the pair production functions
    to set the component of the primary;
    sets the parameters \f$v_{Min}=4 \frac{m_e}{E_p}\f$,
    \f$v_{Max}= max(-\frac{3}{4}\sqrt{e}\frac{m_p}{E_p}Z^{\frac{1}{3}},
    1-6\Big(\frac{m_p}{E_P}\Big)^2, 1-\frac{m_p}{E_p} )\f$
    and \f$v_{Up}=max( v_{max} , v_{Cut}(E_p))\f$
    */

    void setEnergy(int i);

//----------------------------------------------------------------------------//
    /*!
    this is the calculation of the d2Sigma/dvdRo - interface to Integral,
    the function which is defined here is:
    \f[ f(r) =return= \alpha^2r_e^2 \frac{2Z}{1,5\pi}(Z+k)
    \frac{1-v}{v}lpm(r^2,b,s)(F_e+\frac{m_e^2}{m_p^2}F_m)\f]
    */

    double function(double r);

//----------------------------------------------------------------------------//


    /*!
    Landau Pomeranchuk Migdal effect evaluation,
    if the Landau Pomeranchuk Migdal effect is considered in
    the calculation, function is modified by a factor
    \f[lpm=return=\frac{(1+b)(A+(1+r^2)B)+b(C+(1+r^2)D)+(1-r^2)E}{[(2+r^2)(1+b)
    +x(3+r^2)]\ln\Big(1+\frac{1}{x}\Big)+\frac{1-r^2-b}{1+x}-(3+r^2)}\f]
    */

    double lpm(double r2, double b, double x);

//----------------------------------------------------------------------------//
    /*!
    this is the calculation of the dSigma/dv:
    \f[e_{Pair}=return = \rho N_Z z^2 \Big[ \int_{1-r_{max}}^{aux}
    f(r)dr + \int^1_{aux} f(r) dr \Big]\f]
    with \f$ aux=max(1-r_{Max}, ComputerPrecision)\f$ and \f$r_{max} =
    \sqrt{1-\frac{4m_e}{E_p v}}\Big(1-\frac{6m_p^2}{E_p^2(1-v)}\Big)\f$
    */

    double ePair(double v, int i);

//----------------------------------------------------------------------------//
    /*!
    2d parametrization - interface to Interpolate
    if \f$v_{Up}=v_{Max}\f$: \f[return=0;\f]
    else: \f$v=v_{Up}\exp\big(v\log\big(\frac{v_{Max}}{v_{Up}}\big)\big)\f$
    \f[return=e_{pair}(v, Num_{comp})\f]

    */
//----------------------------------------------------------------------------//
    using CrossSections::functionInt;
    double functionInt(double e, double v);

//----------------------------------------------------------------------------//

    // Getter

    double get_vMin()
    {
        return vMin;
    }

    double get_vUp()
    {
        return vUp;
    }

    EpairContinuous* get_Continuous()
    {
        return continuous_;
    }

    EpairStochastic* get_Stochastic()
    {
        return stochastic_;
    }
    
    virtual void SetRandomNumberGenerator(boost::function<double ()> &);

//----------------------------------------------------------------------------//
    //Setter

    void set_jt(bool newJT)
    {
        jt_ =   newJT;
    }
};


#endif /* EPAIRPRODUCTION_H_ */
