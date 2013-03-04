/*! \file   PhotoStochastic.h
*   \brief  Header file for the stochastic photonuclear routines.
*
*   For more details see the class documentation.
*
*   \date   23.06.2010
*   \author Jan-Hendrik Koehne
*
*/

#ifndef PHOTOSTOCHASTIC_H_
#define PHOTOSTOCHASTIC_H_

#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Interpolate.h"
#include "PROPOSAL/Photonuclear.h"


/*! \class PhotoStochastic PhotoStochastic.h "PhotoStochastic.h"
   \brief class contains functions necessary for calculation of photonuclear losses
 */


class PhotoStochastic: public Photonuclear
{

protected:

    std::vector<Integral *> integral_;
    double  *H_;
    double  sum_;
    bool    jt_;
    double  rnd;



//----------------------------------------------------------------------------//

public:

    Interpolate *interpolateJ_;
    Interpolate *interpolateJo_;

//----------------------------------------------------------------------------//
    /*!
     * \brief initializes an integral class separate from that in Propagate
     *
     * \param       *cros   photonucealr crossection
     */

    PhotoStochastic(Photonuclear *cros);

//----------------------------------------------------------------------------//
    /*!
     * this is what d2N/dv,dx is equal to - interface to Integral \n
     * \f$\frac{d^2N}{dv\,dx} = c_p\cdot \frac{d\sigma}{dv} \f$
     *
     * \param   v   relative energy loss
     * \return  \f$\frac{d^2N}{dv\,dx}\f$

     */

    double function(double v);
//----------------------------------------------------------------------------//

    /*!
     * this is what the sum over all components of dN/dx is equal to\n
     * \f[ \frac{dN}{dx} = c_p E_p \sum_{numC}
     * \int_{v_{Min}}^{v_{Up}} v \frac{d\sigma}{dv} dv \f]
     *
     * \return \f$\frac{dN}{dx}\f$
     */

    double dNdx();

//----------------------------------------------------------------------------//
    /*!
     * this is what the sum over all components of dN/dx is equal to
     *
     * \param   rnd     random number which the class member rnd will be set
     * \return  sum over all components
     */

    double dNdx(double rnd);

//----------------------------------------------------------------------------//
    /*!
     * this is the value of \f$e=v \cdot E\f$,
     * corresponding to rnd in the call to dNdx
     * for the component chosen with rnd in the argument
     *
     * \param   rnd     random number to determine on which
     *                  nucleon will be scattered
     * \return  energy loss
     */

    double e(double rnd);

//----------------------------------------------------------------------------//
    /*!
     * 2d parametrization - interface to Interpolate;
     * if \f$v_{Up}=v_{Max}\f$, return 0;
     * else: \f$v=v_{Up}\exp\big(v\log\big(\frac{v_{Max}}{v_{Up}}\big)\big)\f$
     * \f[return=\int_{v_{up}}^{v}photoN dv\f]
     *
     * \param   e   particle will be set to energy e
     * \param   v   relative energy loss
     * \return  0   OR  \f$\int_{v_{up}}^{v}photoN dv\f$
    */

    double functionInt(double e, double v);
//----------------------------------------------------------------------------//
    /*!
     * 1d parametrization - interface to Interpolate
     *
     * \param   e   particle energy e
     * \return  interpolated range calculation for active component
     */

    double functionInt(double e);

//----------------------------------------------------------------------------//
    //Setter

    void set_jt(bool newJT)
    {
        jt_ = newJT;
    }


};
#endif /* PHOTOSTOCHASTIC_H_ */
