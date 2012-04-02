/*
 * Bremsstrahlung.h
 *
 *  Created on: 21.06.2010
 *      Author: koehne
 */

#ifndef Bremsstrahlung_H
#define Bremsstrahlung_H


#include "Integral.h"
#include <algorithm>
#include <cmath>
#include "CrossSections.h"


/*! \class Bremsstrahlung Bremsstrahlung.h "Bremsstrahlung.h"
    \brief class contains functions necessary for calculation of bremsstrahlung losses
 */

class BremsContinuous;

class BremsStochastic;



class Bremsstrahlung: public CrossSections {


        protected:


		double vMax_;			// Set to 0 in Constructor
		double vUp_;			// Set to 0 in Constructor
		double vMin_;			// Set to 0 in Constructor

		int form_;			// Set to 1 in Constrcutor - Switches the parametrization which is used.

		bool init_;		// set to true in Constructor
		double eLpm_;
		double xo_;

		bool lorenz_;		// Set to false in constructor		
		double lorenzCut_;  	/// in [MeV] // - set to 1.e6 in Constructor

                Bremsstrahlung *brems_;

                BremsContinuous *continuous_;
                BremsStochastic *stochastic_;

                // The following functions define the different parametrization for the Bremsstrahlung Cross Section

                double Kelner_Kakoulin_Petrukhin_parametrization(double v, int i);
                double Andreev_Bezrukov_Bugaev_parametrization(double v, int i);
                double Petrukhin_Shestakov_parametrization(double v, int i);
                double complete_screening_case(double v, int i);

	public:



		// constructors

                /*!
                Create an Bremsstrahlung.
                */

                Bremsstrahlung();

                /*!
                creates internal references to p and m, to be called from subclasses
                */

                Bremsstrahlung(Bremsstrahlung *cros);

                /*!
                initializes subclasses and creates internal references to p and m
                */

                Bremsstrahlung(CrossSections *cros);


		// Memberfunctions
                //---------------------------------------------------------------------------------------------------------//

                /*!
                call before using the bremsstrahlung functions to set the component of the primary;
                \f[v_{Max}=1-\frac{3}{4}\sqrt{e}\frac{m}{E}Z^{\frac{1}{3}}\f];
                \f[ v_{Max}=min(v_{Max}, 1-\frac{m}{E})\f];
                \f[v_{up}=min(v_{Max}, v_{Cut})\f]
                */

		void setEnergy(int i);

                //---------------------------------------------------------------------------------------------------------//

                /*!
                this is what the Elastic Bremsstrahlung Cross Section (EBCS) is equal to
                units are [1/cm] since the multiplication by No*n is done here.
                Corrections for excitations of the nucleus and deep inelastic excitations of separate nucleons are
                included (positive term dn/Z), as well as the contribution of the mu-diagrams to the inelastic
                bremsstrahlung on the electrons (non-zero only for allowed energies of photon after electron recoil).
                four different parametrizations are enumerated, the final result is:
                \f[ \sigma_{el}= \rho_{mol}N(z^2)^2\Big[a_1\Big] \f] and $a_1$ depends on the chosen parametrization
                */

                double Sel(double v, int i);

                //---------------------------------------------------------------------------------------------------------//

                /*!
                Landau Pomeranchuk Migdal effect and dielectric suppression evaluation
                \f$ e_{LPM} \f$ is evaluated:
                \f[e_{LPM}=\frac{(\alpha m)^2}{4\pi m_e r_e \sum()} \f]
                with \f$\sum()=\sum_{numComp}\Big(\int_0^{v_{up}}v\sigma dv +\int_{v_{up}}^{v_{Max}}v\sigma dv \Big)\f$
                */

                void setLpm();

                //---------------------------------------------------------------------------------------------------------//

                /*!
                Landau Pomeranchuk Migdal effect and dielectric suppression evaluation
                lpm is evaluated:
                \f[ lpm= return=\frac{x_i}{3}\Big[v^2\frac{G(s)}{\gamma^2}+2(1+(1-v)^2)\frac{\Phi(s)}{\gamma}\Big] \f]
                */

                double lpm(double v, double s1);

                //---------------------------------------------------------------------------------------------------------//


                // Getter

                        BremsContinuous* get_Continuous();
                        BremsStochastic* get_Stochastic();
                        Bremsstrahlung* get_Bremsstrahlung() {return brems_;}

			double get_vMax() const {return vMax_;}
			double get_vUp() const {return vUp_;}
			double get_vMin() const {return vMin_;}

                        Integral* get_integral() const {return integral_;}

			int get_form() const {return form_;}

			bool get_init() const {return init_;}
			double get_eLpm() const {return eLpm_;}
			double get_xo() const {return xo_;}

                        //bool get_lorenz() const {return lorenz_;};
                        //double get_lorenzCut() const {return lorenzCut_;}

			// Setter

			void set_continous(BremsContinuous *continous);
			void set_stochastic(BremsStochastic *stochastic);

			void set_vMax(double vMax);
			void set_vUp(double vUp);
			void set_vMin(double vMin);

			void set_integral(Integral *integral);

			void set_form(int form);

			void set_init(bool init);
			void set_eLpm(double eLpm);
			void set_xo(double xo);

			void set_lorenz(bool lorenz);
			void set_lorenzCut(double lorenzCut);


		// destructors

		///@brief Crush this Bremsstrahlung.
		virtual ~Bremsstrahlung();


};



#endif //Bremsstrahlung_H
