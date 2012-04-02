/*
 * IonizContinuous.h
 *
 *  Created on: 24.06.2010
 *      Author: koehne
 */

#ifndef IONIZCONTINUOUS_H_
#define IONIZCONTINUOUS_H_

#include "Integral.h"
#include "Ionizationloss.h"
#include "Interpolate.h"

/*!
  \brief class contains functions necessary for calculation of ionization losses
 */


class IonizContinuous : public Ionizationloss {

	protected:

		Integral *integral_;


                bool jt_;                    // Set to false in conststructor

		//----------------------------------------------------------------------------------------------------//

                /*!
                 initializes an integral class separate from that in Propagate
		 */

	public:


                Interpolate *interpolateJ_;

		IonizContinuous(Ionizationloss *cros);

		//----------------------------------------------------------------------------------------------------//

                /*!
                 inelastic electron bremsstrahlung correction to dEdx - interface to Integral;
                \f[f(r)= v \frac{d^2 N}{dv dx} \cdot \Delta_{inel}\f]
		 */

		double function(double v);
		//----------------------------------------------------------------------------------------------------//

                /*!
                 density correction term calculation;
                 */
                 /*!
                 if \f$X=\frac{\ln\frac{\gamma}{\beta}}{\ln10} < X_0\f$  :  \f$\delta= \delta_0 10^{2(X-X_0)}\f$
                */
                /*!
                if \f$ X<X_1\f$ :  \f$\delta= c_1X+c+a(X1-X)^m\f$
                */
                /*!
                either :  \f$\delta = c_1x+c  \f$
                 */

		double delta();

		//----------------------------------------------------------------------------------------------------//

                /*!
                  contribution of ionization to -dE/dx;
                  \f[\frac{dE'}{dx}=return=c_i\Big[\rho \cdot \frac{dE}{dx} + E \int_{v_{min}}^{v_{up}} v \frac{d^2N}{dvdx}\Delta_{inelastic}dv\Big]\f]
                  with
                  \f[ \frac{dE}{dx} =Kz^2\frac{Z}{A}\frac{1}{2\beta^2}\Big[ \ln(2v_{up}m_e E) +2\ln\Big(\frac{\beta \gamma}{10^{-6}I}\Big)+\Big(\frac{v_{up}}{2\Big(1+\frac{1}{\gamma}\Big)} \Big)^2-\beta^2\Big(1+\frac{v_{up}}{v_{Max}}\Big)-\delta \Big] \f]

                */

		double dEdx();

		//----------------------------------------------------------------------------------------------------//


                /*!
                 1d parametrization - interface to Interpolate
		 */

                double functionInt(double e);

                // Setter
                void set_jt(bool newJT);
};

#endif /* IONIZCONTINUOUS_H_ */
