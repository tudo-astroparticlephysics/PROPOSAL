/*
 * PhotoStochastic.h
 *
 *  Created on: 23.06.2010
 *      Author: koehne
 */

#ifndef PHOTOSTOCHASTIC_H_
#define PHOTOSTOCHASTIC_H_

#include "Integral.h"
#include "Interpolate.h"
#include "Photonuclear.h"


/*! \class PhotoStochastic PhotoStochastic.h "PhotoStochastic.h"
   \brief class contains functions necessary for calculation of photonuclear losses
 */


class PhotoStochastic: public Photonuclear{

	protected:

                std::vector<Integral*> integral_;
                double* H_;
		double sum_;

                bool jt_;		// set to false in constructor
		double rnd;



		//----------------------------------------------------------------------------------------------------//

                /*!
                  initializes an integral class separate from that in Propagate
		 */
	public:

                Interpolate* interpolateJ_;
                Interpolate* interpolateJo_;

		PhotoStochastic(Photonuclear *cros);

		//----------------------------------------------------------------------------------------------------//

                /*!
                  this is what d2N/dv,dx is equal to - interface to Integral;
                  \f$return=c_p\cdot \frac{d\sigma}{dv} \f$
		 */

		double function(double v);
		//----------------------------------------------------------------------------------------------------//

                /*!
                  this is what the sum over all components of dN/dx is equal to;
                  \f[ \frac{dN}{dx} = c_p E_p \sum_{numC} \int_{v_{Min}}^{v_{Up}} v \frac{d\sigma}{dv} dv \f]
		 */

		double dNdx();

		//----------------------------------------------------------------------------------------------------//

                /*!
                  this is what the sum over all components of dN/dx is equal to
		 */

		double dNdx(double rnd);

		//----------------------------------------------------------------------------------------------------//

                /*!
                  this is the value of \f$e=v \cdot E\f$, corresponding to rnd in the call to dNdx
                  for the component chosen with rnd in the argument
		 */

		double e(double rnd);

		//----------------------------------------------------------------------------------------------------//

                /*!
                  2d parametrization - interface to Interpolate;
                  if \f$v_{Up}=v_{Max}\f$, return 0;
                  else: \f$v=v_{Up}\exp\big(v\log\big(\frac{v_{Max}}{v_{Up}}\big)\big)\f$
                  \f[return=\int_{v_{up}}^{v}photoN dv\f]
		 */

		double functionInt(double e, double v);
		//----------------------------------------------------------------------------------------------------//

                /*!
                  1d parametrization - interface to Interpolate
		 */

		double functionInt(double e);


                //Setter

                void set_jt(bool newJT){
                    jt_ = newJT;
                }
};
#endif /* PHOTOSTOCHASTIC_H_ */
