/*
 * Scattering.h
 *
 *  Created on: 29.06.2010
 *      Author: koehne
 */

#ifndef SCATTERING_H_
#define SCATTERING_H_

class Interpolate;
class PROPOSALParticle;
class Propagate;

#include "Integral.h"
#include "PhysicsModel.h"
#include "Output.h"


/*! \class Scattering Scattering.h "Scattering.h"
   \brief class contains functions for evaluation of the Moliere angle distribution spread
 */


class Scattering : public PhysicsModel {

protected:

        PROPOSALParticle *particle_;
        Propagate *propagate_;
        Integral *integral_;
         Output *o;

        static double cutoff;

    	bool df; 		// Preset as false in constructor

	bool jt; 	// Preset as false in constructor


    //----------------------------------------------------------------------------------------------------//

    /**
     * initialize the class
     */
public:
        Interpolate *interpolateJ_;
        Interpolate *interpolateJdf_;

	Scattering(PROPOSALParticle *p);

    //----------------------------------------------------------------------------------------------------//

    /*!
      function for angle distribution spread calculation - interface to Integral;
     \f[f(E) = \Big( R_y \frac{E_p}{p_p^2}\Big)^2 \frac{1}{  \frac{dE}{dx}\big|_{Ioniz} +\frac{dE}{dx}\big|_{Brems}+\frac{dE}{dx}\big|_{PhotoNuc}+\frac{dE}{dx}\big|_{epair}  }\f]
     */

	double function(double E);

    //----------------------------------------------------------------------------------------------------//

    /*!
     spread of the angle distribution, corresponding to the given propagation distance;
     \f[\theta_0 = \frac{\sqrt{\int_{e_i}^{e_f} \Big( R_y \frac{E_p}{p_p^2}\Big)^2 \frac{1}{  \frac{dE}{dx}\big|_{Ioniz} +\frac{dE}{dx}\big|_{Brems}+\frac{dE}{dx}\big|_{PhotoNuc}+\frac{dE}{dx}\big|_{epair}  } }}{X_0} z \Big(1+0,038\ln\big(\frac{dr}{X_0}\big)\Big)\f]
     */

	double gettho(double dr, double ei, double ef);

    //----------------------------------------------------------------------------------------------------//

    /**
     * 1d parametrization - interface to Interpolate
     */

	double functionInt(double e);




// Getter

    bool get_jt(){return jt;}

    bool get_df(){return df;}

    Interpolate* get_interpolateJ(){return interpolateJ_;}

    Interpolate* get_interpolateJdf(){return interpolateJdf_;}
    Propagate* get_propagate(){return propagate_;}

// Setter

    void set_jt(bool newJT );

    void set_df(bool newDF );
};

#endif /* SCATTERING_H_ */
