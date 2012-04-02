/*
 * Ionizationloss.h
 *
 *  Created on: 22.06.2010
 *      Author: koehne
 */

#ifndef IONIZATIONLOSS_H_
#define IONIZATIONLOSS_H_

#include "CrossSections.h"

/*!
\class Ionizationloss Ionizationloss.h "Ionizationloss.h"
\brief class contains functions necessary for calculation of ionization losses
*/



class IonizContinuous;

class IonizStochastic;

class Ionizationloss : public CrossSections{

	protected:


		double vMax, vUp, vMin;
		double beta, gamma;
		Ionizationloss *ioniz;
                IonizContinuous *cont;
                IonizStochastic *stochastic_;


	public:


	// Called by Subclasses

        //Ionizationloss();

    /*!
      creates internal references to p and m, to be called from subclasses
     */

	Ionizationloss(Ionizationloss *cros);

    /*!
      initializes subclasses and creates internal references to p and m
     */

	Ionizationloss(CrossSections *cros);

//---------------------------------------------------------------------------------------------------------//

    /*!
      call before using the ionization functions, but after setting the particle energy;
      sets the parameters \f$v_{min}=\frac{10^{-6}I}{E_p}\f$, \f$v_{Max}=min\Big(v_{Max}, 1-\frac{m_p}{E_p}\Big)\f$
      with \f$v_{Max}=\frac{2m_e(\gamma^2-1)}{1+2\gamma\frac{m_e}{m_p}+\Big(\frac{m_e}{m_p}\Big)^2E}\f$ and
      \f$v_{up}=min(v_{Max}, v_{Cut})\f$
     */

	void setEnergy();

//---------------------------------------------------------------------------------------------------------//

    /*!
      this is what d2Ndvdx is equal to in the first approximation:
      \f[\frac{d^2N}{dxdv}=return =\frac{Kz^2Z\rho}{2A\beta^2Ev^2} \Big[1-\beta^2 \frac{v}{v_{max}}+\frac{1}{2}\Big(\frac{v}{1+\frac{1}{\gamma}}\Big)^2 \Big]\f]
    */

	double d2Ndvdx(double v);

//---------------------------------------------------------------------------------------------------------//

    /*!
      this is the inelastic electron bremsstrahlung correction to the first approximation of the \f$ \frac{d2N}{dvdx}\f$:
      \f[\Delta_{inelastic}=return =\frac{\alpha}{2\pi} \Big[ \ln\big(1+\frac{2vE}{m_e}\big)\Big[ 2\ln\Big(\frac{1-\frac{v}{v_{max}}}{1-v}\Big)+\ln\Big(\frac{2\gamma m_e(1-v)}{mv}\Big)\Big]-\ln^2\Big(\frac{1-\frac{v}{v_{max}}}{1-v}\Big) \Big]\f]
     */

	double inelCorrection(double v);

//---------------------------------------------------------------------------------------------------------//

	//Getter
        IonizContinuous* get_Continuous();
        IonizStochastic* get_Stochastic();
	double get_beta() { return beta; }
        double get_vMin() {return vMin;}
        double get_vMax() {return vMax;}
        double get_vUp()  {return vUp;}
};
#endif /* IONIZATIONLOSS_H_ */
