/*
 * IonizStochastic.h
 *
 *  Created on: 24.06.2010
 *      Author: koehne
 */

#ifndef IONIZSTOCHASTIC_H_
#define IONIZSTOCHASTIC_H_


#include "Ionizationloss.h"


/*!
 \brief class contains functions necessary for calculation of ionization losses
 */

class IonizStochastic: public  Ionizationloss {

	protected:
		Integral *integral_;

                bool jt_; 			// Set to false in constructor
		double rnd, sum;



		//----------------------------------------------------------------------------------------------------//

                /*!
                  initializes an integral class separate from that in Propagate
		 */

	public:
                Interpolate *interpolateJ_;
                Interpolate *interpolateJo_;



		IonizStochastic(Ionizationloss *cros);

		//----------------------------------------------------------------------------------------------------//

                /*!
                 this is what d2Ndvdx is equal to - interface to Integral;
                  \f[ return = c_i \frac{d^2N}{dvdx} (1+\Delta_{inel}) \f]
                 */

		double function(double v);

		//----------------------------------------------------------------------------------------------------//

                /*!
                 this is what dNdx is equal to:
                 \f[\frac{dN}{dx} =return=c_i \cdot \int \frac{d^2N}{dvdx} (1+\Delta_{inelastic})dv\f]
		 */

		double dNdx();

		//----------------------------------------------------------------------------------------------------//

                /*!
                 this is what dNdx is equal to;
                 like dNdx(), but a double rnd is assigned, which is important for getUpper
		 */

		double dNdx(double rnd);

		//----------------------------------------------------------------------------------------------------//

                /*!
                  this is the value of \f$e=v*E_{particle}\f$, corresponding to rnd in the call to dNdx
		 */

		double e(double rnd);

		//----------------------------------------------------------------------------------------------------//


                /*!
                 2d parametrization - interface to Interpolate
                 if \f$v_{Up}=v_{Max}\f$: \f[return=0;\f]
                 else: \f$v=v_{Up}\exp\big(v\log\big(\frac{v_{Max}}{v_{Up}}\big)\big)\f$
                 \f[return=\int_{v_{up}}^{v} c_i \frac{d^2N}{dvdx} (1+\Delta_{inel}) dv\f]

		 */

		double functionInt(double e, double v);

		//----------------------------------------------------------------------------------------------------//



                /*!
                 1d parametrization - interface to Interpolate
		 */

		double functionInt(double e);


                // Setter

                void set_jt(bool newJT);

                // Getter

                Interpolate* get_interpolateJ();
                Interpolate* get_interpolateJo();
};


#endif /* IONIZSTOCHASTIC_H_ */
