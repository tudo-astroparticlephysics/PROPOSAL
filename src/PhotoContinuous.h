/*
 * PhotoContinuous.h
 *
 *  Created on: 28.06.2010
 *      Author: koehne
 */

#ifndef PHOTOCONTINUOUS_H_
#define PHOTOCONTINUOUS_H_

#include "Photonuclear.h"
#include "Integral.h"
#include "Interpolate.h"



/*! \class PhotoContinuous PhotoContinuous.h "PhotoContinuous.h"
   \brief class contains functions necessary for calculation of photonuclear losses
 */




class PhotoContinuous: public Photonuclear{

protected:
	Integral *integral_;

        bool jt_;		// Set to false in constructor

    //----------------------------------------------------------------------------------------------------//

    /*!
      initializes an integral class separate from that in Propagate
     */

public:

        Interpolate *interpolateJ_;


	PhotoContinuous(Photonuclear *cros);

    //----------------------------------------------------------------------------------------------------//

    /*!
      photonuclear energy losses - interface to Integral;
      \f[return= v\cdot photoN\f]
     */

	double function(double v);

    //----------------------------------------------------------------------------------------------------//

    /*!
      Contribution of photonuclear interactions to -dE/dx;
      \f[\frac{dE}{dx}=return = c_p E_p \sum_{numC} \int_{v_{Min}}^{v_{Up}} v \frac{d\sigma}{dv} dv\f]
     */

	double dEdx();

    //----------------------------------------------------------------------------------------------------//

    /*!
      1d parametrization - interface to Interpolate
     */

	double functionInt(double e);

        void set_jt(bool newJT)
        {
            jt_ = newJT;
        }

   // Destructor

        ~PhotoContinuous(){};

};


#endif /* PHOTOCONTINUOUS_H_ */
