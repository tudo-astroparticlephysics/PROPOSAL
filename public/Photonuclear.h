/*
 * Photonuclear.h
 *
 *  Created on: 29.06.2010
 *      Author: koehne
 */

#ifndef PHOTONUCLEAR_H_
#define PHOTONUCLEAR_H_

#include "CrossSections.h"

#include "Integral.h"
#include "Interpolate.h"
#include <vector>


/*! \class Photonuclear Photonuclear.h "Photonuclear.h"
   \brief class contains functions necessary for calculation of photonuclear losses
 */



class PhotoContinuous;
class PhotoStochastic;


class Photonuclear: public  CrossSections{

	protected:
                PhotoContinuous *continuous_;
                PhotoStochastic *stochastic_;

		double vMax, vUp, vMin;
		Photonuclear *photo;

		Integral *integral_;
		int i;
		double v;
		int bb;				// set to 1 in constructor
                int form;			// set to 1 in constructor

		bool initM;			// set to true in Constructor
		Interpolate *interpolateM_;

		bool initH;			// Set to true in Constructor
		const static int hmax=8;
                Interpolate* interpolateH_;
		int shadow;			// set to 1 in constructor

                bool jt_;			// set to false in constructor

                //----------------------------------------------------------------------------------------------------//

                /*!
                  creates internal references to p and m, to be called from subclasses
		 */

	public:

                Interpolate* interpolateJ_;
		Photonuclear(Photonuclear *cros);

		//----------------------------------------------------------------------------------------------------//

                /*!
                  initializes subclasses and creates internal references to p and m
		 */

		Photonuclear(CrossSections *cros);

		//----------------------------------------------------------------------------------------------------//

                /*!
                   call before using the photonuclear functions to set the component of the primary
                   sets the parameters \f$v_{min}=\frac{1}{E_{p}}\Big(m_{\pi}+\frac{m_{\pi}^2}{2M_N}\Big)\f$, \f$v_{Max}=min\Big(v_{Max}, 1-\frac{m_p}{E_p}\Big)\f$
                   with \f$v_{Max}=1-\frac{M_N\Big[1+\Big(\frac{m_p}{M_N}\Big)^2\Big]}{2E_p}\f$ if \f$m_p<m_{\pi}\f$ otherwise
                   \f$ v_{Max}=1\f$;
                   \f$v_{up}=min(v_{Max}, v_{Cut}(E))\f$
                */

		void setEnergy(int i);

		//----------------------------------------------------------------------------------------------------//

                /*!
                  this is what the photonuclear interaction cross section is equal to,
                  there are different cases in which the cross section \f$ \sigma_{\gamma,N}\f$ is calculated,
                  intermediate result is:
                  \f[f(v) = \frac{\alpha}{2\pi}A \sigma_{\gamma, N}v\Big[\frac{3}{4}G \Big[k\ln\Big(1+\frac{m_1}{t}\Big)+ \frac{4\mu^2}{m_1}\ln\Big(1+\frac{m_1}{t}\Big)-\frac{km_1}{m_1+t}-\frac{2\mu^2}{t} \Big] + \frac{1}{4}\Big[\Big(k+\frac{2\mu^2}{m_2}\Big)\ln\Big(1+\frac{m_2}{t}\Big)-\frac{2\mu^2}{t}\Big]
                  + \frac{1}{2}\frac{\mu^2}{t}\Big[\frac{3}{4}G\frac{m_1-4t}{m_1+t}+\frac{1}{4}\frac{m_2}{t}\ln\Big(1+\frac{t}{m_2}\Big)\Big] \Big]\f]
                  the final result is:
                  \f[return=\rho_{mol}N(z^2)^2 \int_{min}^{max} \frac{d\sigma}{dv} dv\f]
                  with
                  \f[\frac{d\sigma}{dv} = \rho_{mol}N(z^2)^2\cdot f(v)\f]
                */

		double photoN(double v, int i);


		//----------------------------------------------------------------------------------------------------//

                /*!
                  call this to replace photon-nucleon cross section approximation formula with
                  a more exact at low energies parametrization of the experimental data
		 */

		void setMeasured();

		//----------------------------------------------------------------------------------------------------//

                /*!
                  returns measured value of the photon nucleon cross section
		 */

		double measuredSgN(double e);

		//----------------------------------------------------------------------------------------------------//

                /*!
                  initializes hard part of BB (by Bugaev and Shlepin) parametrization calculation
		 */

		void enableHardBB();

		//----------------------------------------------------------------------------------------------------//

                /*!
                  returns interpolated value of the hard part of bb cross section
		 */

		double hardBB(double e, double v);

		//----------------------------------------------------------------------------------------------------//

                /*!
                  parametrized photonuclear cross section - interface to Integral
		 */

		double function(double Q2);

		//----------------------------------------------------------------------------------------------------//

                /*!
                  2d parametrization - interface to Interpolate;
                  if \f$v_{up}=v_{Max}\f$ return=0; otherwise set
                  \f[v=v_{Up}\cdot \exp \big(v \ln \big[\frac{v_{Max}}{v_{up}}\big]\big)\f]
                  and return the result from Photonuclear::photoN with different v
		 */

		double functionInt(double e, double v);

		// Getter

                PhotoContinuous* get_Continuous();
                PhotoStochastic* get_Stochastic();
                double get_vMin(){return vMin;}
                double get_vUp(){return vUp;}
                double get_vMax(){return vMax;}
                int get_form() {return form;}
                bool get_jt() {return jt_;}
                int get_bb(){return bb;}
                int get_shadow(){return shadow;}


                // Setter

               void set_jt(bool newjt);
               void set_form(int newForm);
               void set_bb(int newbb);
               void set_shadow(int newshadow);



};


#endif /* PHOTONUCLEAR_H_ */
