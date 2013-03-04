/*! \file   IonizStochastic.h
*   \brief  Header file for the stochastic ionizationloss routines.
*
*   For more details see the class documentation.
*
*   \date   24.06.2010
*   \author Jan-Hendrik Koehne
*/

#ifndef IONIZSTOCHASTIC_H_
#define IONIZSTOCHASTIC_H_


#include "PROPOSAL/Ionizationloss.h"


/*!
 *\brief class contains functions necessary for calculation of ionization losses
 */

class IonizStochastic: public  Ionizationloss
{

protected:

    Integral    *integral_;         //!< integral class
    bool        jt_;                //!< interpolate results?;
    double      rnd;
    double      sum;



//----------------------------------------------------------------------------//


public:

    Interpolate *interpolateJ_;     //!< interpolate object for integrated results
    Interpolate *interpolateJo_;    //!< interpolate object for differential function

//----------------------------------------------------------------------------//
    /*!
     * initializes an integral class separate from that in Propagate
     *
     * \param   *cros   Ionization crossection for initializing stochastic losses
     */

    IonizStochastic(Ionizationloss *cros);

//----------------------------------------------------------------------------//
    /*!
     * this is what d2Ndvdx is equal to - interface to Integral.
     *
     * \return \f[ return = c_i \frac{d^2N}{dvdx} (1+\Delta_{inel}) \f]
     */

    double function(double v);

//----------------------------------------------------------------------------//
    /*!
     * this is what dNdx is equal to:
     * \f[\frac{dN}{dx} =c_i \cdot \int
     * \frac{d^2N}{dvdx} (1+\Delta_{inelastic})dv\f]
     *
     * \return dN/dx
     */

    double dNdx();

//----------------------------------------------------------------------------//
    /*!
     * this is what dNdx is equal to;
     * like dNdx(), but a double rnd is assigned,
     * which is important for getUpper
     *
     * \param   rnd     important for getUpper
     * \return  dN/dx
     */

    double dNdx(double rnd);

//----------------------------------------------------------------------------//
    /*!
     * this is the value of \f$e=v*E_{particle}\f$,
     * corresponding to rnd in the call to dNdx
     *
     * \param   rnd
     * \return  energyloss by inionization
     */

    double e(double rnd);

//----------------------------------------------------------------------------//

    /*!
     * 2d parametrization - interface to Interpolate
     * if \f$v_{Up}=v_{Max}\f$: \f[return=0;\f]
     * else: \f$v=v_{Up}\exp\big(v\log\big(\frac{v_{Max}}{v_{Up}}\big)\big)\f$
     * \f[return=\int_{v_{up}}^{v} c_i \frac{d^2N}{dvdx} (1+\Delta_{inel}) dv\f]
     *
     * \param   e   particle energy
     * \param   v   relative energy loss
     * \return
    */

    double functionInt(double e, double v);

//----------------------------------------------------------------------------//


    /*!
     * 1d parametrization - interface to Interpolate
     *
     * \param   v   relative energy loss
     * \return
     */

    double functionInt(double e);

//----------------------------------------------------------------------------//
    // Setter

    void set_jt(bool newJT);

//----------------------------------------------------------------------------//
    // Getter

    Interpolate* get_interpolateJ();
    Interpolate* get_interpolateJo();
};


#endif /* IONIZSTOCHASTIC_H_ */
