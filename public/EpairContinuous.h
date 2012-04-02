/*
 * EpairContinuous.h
 *
 *  Created on: 22.06.2010
 *      Author: koehne
 */

#ifndef EPAIRCONTINUOUS_H_
#define EPAIRCONTINUOUS_H_

#include "Epairproduction.h"
#include "Integral.h"
#include "Interpolate.h"



/*! \class EpairContinuous EpairContinuous.h "EpairContinuous.h"
   \brief class contains functions necessary for calculation of e+e- pair production losses
 */


class EpairContinuous: public Epairproduction{

        protected:
                Integral *integral_;
		bool reverse_;			// Set to false in Constructor

                bool jt_;			// Set to false in Constructor

    //----------------------------------------------------------------------------------------------------//


	public:


                Interpolate *interpolateJ_;

	// Default Constructor


		EpairContinuous() ;

    /*!
     * initializes an integral class separate from that in Propagate
     */

		EpairContinuous(Epairproduction *cros);


    /*!
      pair production energy losses - interface to Integral;
      \f[ return= v\cdot e_{Pair}(v, components)\f]
     */

		double function(double v);

    //----------------------------------------------------------------------------------------------------//

    /*!
     contribution of pair production to -dE/dx,
     \f[\frac{dE}{dx}=return= c_e E_p \sum_{numC} \int_{v_{Min}}^{v_{Up}} v e_{Pair}(v, components) dv\f]
    */

		double dEdx();

    //----------------------------------------------------------------------------------------------------//

    /*!
      1d parametrization - interface to Interpolate
     */

		double functionInt(double e);

    // Getter

                Interpolate* get_J(){ return interpolateJ_;}

                void set_jt(bool newjt){jt_=newjt;}

};


#endif /* EPAIRCONTINUOUS_H_ */
