/*! \file   Scattering.h
*   \brief  Header file for the scattering routines.
*
*   For more details see the class documentation.
*
*   \date   29.06.2010
*   \author Jan-Hendrik Koehne
*/

#ifndef SCATTERING_H_
#define SCATTERING_H_

class Interpolate;
class PROPOSALParticle;
class Propagate;

#include "PROPOSAL/Integral.h"
#include "PROPOSAL/PhysicsModel.h"
#include "PROPOSAL/Output.h"


/*! \class Scattering Scattering.h "Scattering.h"
   \brief class contains functions for evaluation
    of the Moliere angle distribution spread
 */


class Scattering : public PhysicsModel
{

protected:

    PROPOSALParticle *particle_;
    Propagate *propagate_;
    Integral *integral_;
    Output *o;

    static double cutoff;

    bool df;
    bool jt;


 //----------------------------------------------------------------------------//

public:

    Interpolate *interpolateJ_;
    Interpolate *interpolateJdf_;

//----------------------------------------------------------------------------//
    /**
     * initialize the class
     */

    Scattering(PROPOSALParticle *p);

//----------------------------------------------------------------------------//
    /*!
     * function for angle distribution spread calculation
     * - interface to Integral;
     * \f[f(E) = \Big( R_y \frac{E_p}{p_p^2}\Big)^2
     * \frac{1}{  \frac{dE}{dx}\big|_{Ioniz} +\frac{dE}{dx}\big|_{Brems}
     * +\frac{dE}{dx}\big|_{PhotoNuc}+\frac{dE}{dx}\big|_{epair}  }\f]
     *
     * \param   E   energy
     * \return  \f$f(E)\f$
     */

    double function(double E);

//----------------------------------------------------------------------------//
    /*!
     * spread of the angle distribution,
     * corresponding to the given propagation distance;
     * \f[\theta_0 = \frac{\sqrt{\int_{e_i}^{e_f} \Big( R_y
     * \frac{E_p}{p_p^2}\Big)^2 \frac{1}{  \frac{dE}{dx}\big|_{Ioniz}
     * +\frac{dE}{dx}\big|_{Brems}+\frac{dE}{dx}\big|_{PhotoNuc}+
     * \frac{dE}{dx}\big|_{epair}  } }}{X_0} z
     * \Big(1+0,038\ln\big(\frac{dr}{X_0}\big)\Big)\f]
     *
     * \param   dr  displacement
     * \param   ei  initial energy
     * \param   ef  final energy
     * \return  \f$\theta_0\f$
     */

    double gettho(double dr, double ei, double ef);

//----------------------------------------------------------------------------//
    /**
     * 1d parametrization - interface to Interpolate
     *
     * \param   E   energy
     * \return  differential or integrated value of function(x)
     */

    double functionInt(double e);


//----------------------------------------------------------------------------//
// Getter

    bool get_jt()
    {
        return jt;
    }

    bool get_df()
    {
        return df;
    }

    Interpolate *get_interpolateJ()
    {
        return interpolateJ_;
    }

    Interpolate *get_interpolateJdf()
    {
        return interpolateJdf_;
    }
    Propagate *get_propagate()
    {
        return propagate_;
    }
//----------------------------------------------------------------------------//
// Setter

    void set_jt(bool newJT );

    void set_df(bool newDF );

};

#endif /* SCATTERING_H_ */
