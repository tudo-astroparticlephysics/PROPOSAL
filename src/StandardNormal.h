/*
 * StandardNormal.h
 *
 *  Created on: 02.08.2010
 *      Author: koehne
 */

#ifndef STANDARDNORMAL_H_
#define STANDARDNORMAL_H_

#include "Integral.h"
#include <cmath>

#include "MathModel.h"
#include "FunctionOfx.h"
#include "FunctionInt.h"
#include "PhysicsModel.h"

class Interpolate;


/**
 * This class provides routines for evaluation of random numbers distributed normally.
 * @author Dmitry Chirkin
 */

class StandardNormal : public MathModel, public FunctionOfx, public FunctionInt{ 

protected:
        PhysicsModel *Physic;
        Integral *integral_;
        double val1, val2;

	bool jt; 				// Setted to false in Constructor
    	double norm;		// = sqrt(2*PI_) in Constructor

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes class with default settings
     */

public:
        Interpolate *interpolateJ_;

	StandardNormal();

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes the class
     */

	StandardNormal(int romberg, int maxSteps, double precision);
    //----------------------------------------------------------------------------------------------------//

    /**
     * class initializer
     */

        void init(int romberg, int maxSteps, double precision);

    //----------------------------------------------------------------------------------------------------//

    /**
     * evaluates the integrated probability
     */

	double sndpr(double x);

    //----------------------------------------------------------------------------------------------------//

    /**
     * evaluates the standard normal random number
     */

	double sndrn(double x);
    //----------------------------------------------------------------------------------------------------//

    /**
     * evaluates the standard normal random number
     */

	double sndrn(double rnd, double average, double sigma, double xmin, double xmax, bool cutoff);

    //----------------------------------------------------------------------------------------------------//

    /**
     * function describes standard normal distribution - interface to Integral
     */

	double function(double x);

    //----------------------------------------------------------------------------------------------------//

    /**
     * 1d parametrization - interface to Interpolate
     */

	double functionInt(double x);


 // Getter

        bool get_jt(){return jt;}

        Interpolate get_interpolateJ();

// Setter

        void set_jt(bool newJT);

};


#endif /* STANDARDNORMAL_H_ */
