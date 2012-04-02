/*
 * Decay.h
 *
 *  Created on: 24.06.2010
 *      Author: koehne
 */

#ifndef DECAY_H_
#define DECAY_H_

#include "CrossSections.h"
#include "FindRoot.h"
#include <string>
#include "Output.h"


/*! \class  Decay  Decay.h " Decay.h"
   \brief class contains functions necessary for the calculation of decay
 */



class Decay: public CrossSections , public DFunctionOfx { 

	protected:

		FindRoot *f;
		bool flag;
		std::string out;


		/**
		 * creates internal references to p and m
		 */

	public:

		Decay(CrossSections *cros);
		//----------------------------------------------------------------------------------------------------//

		/**
		 * this cross section describes decay
		 */

		double decay();

		//----------------------------------------------------------------------------------------------------//


		/**
		 * energy of the electron that results from the muon decay
		 */

		double e(double ernd, double arnd, double srnd, Output *o);

		//----------------------------------------------------------------------------------------------------//

		/**
		 * function for electron energy calculation - interface to FindRoot
		 */

		double function(double x);

		//----------------------------------------------------------------------------------------------------//

		/**
		 * function for electron energy calculation - interface to FindRoot
		 */

		double dFunction(double x);

		// GETTER


		std::string get_out(){ return out; }

};


#endif /* DECAY_H_ */
