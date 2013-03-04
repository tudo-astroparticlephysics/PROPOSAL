/*! \file   StandardNormal.h
*   \brief  Header file for the standard normal routines.
*
*   For more details see the class documentation.
*
*   \date   02.08.2010
*   \author Jan-Hendrik Koehne
*/


#ifndef STANDARDNORMAL_H_
#define STANDARDNORMAL_H_

#include "PROPOSAL/Integral.h"
#include <cmath>

#include "PROPOSAL/MathModel.h"
#include "PROPOSAL/FunctionOfx.h"
#include "PROPOSAL/FunctionInt.h"
#include "PROPOSAL/PhysicsModel.h"

class Interpolate;


/**
 * \brief This class provides routines for evaluation of
 * random numbers distributed normally.
 */

class StandardNormal : public MathModel, public FunctionOfx, public FunctionInt
{

protected:
    PhysicsModel *Physic;
    Integral *integral_;
    double val1, val2;

    bool jt;
    double norm;        	// = sqrt(2*PI_) in Constructor

//----------------------------------------------------------------------------//


public:

    Interpolate *interpolateJ_;

//----------------------------------------------------------------------------//
    //Constructors

    /**
     * \brief initializes class with default settings
     */

    StandardNormal();

//----------------------------------------------------------------------------//
    /**
     * \brief initializes the class
     *
     * \param   romberg     order of interpolation and integration
     * \param   maxSteps    number of sampling points
     * \param   precision   integration precision
     */

    StandardNormal(int romberg, int maxSteps, double precision);
//----------------------------------------------------------------------------//
    //Memberfunctions
    /**
     * \brief class initializer
     *
     * \param   romberg     order of interpolation and integration
     * \param   maxSteps    number of sampling points
     * \param   precision   integration precision
     */

    void init(int romberg, int maxSteps, double precision);

//----------------------------------------------------------------------------//
    /**
     * \breif evaluates the integrated probability
     *
     * \param   x   wanted value
     * \return  integrated probability
     */

    double sndpr(double x);

//----------------------------------------------------------------------------//
    /**
     * evaluates the standard normal random number
     *
     * \param   x   wanted value
     * \return  standard normal random number
     */

    double sndrn(double x);
//----------------------------------------------------------------------------//
    /**
     * evaluates the standard normal random number
     *
     * \param   rnd         random number
     * \param   average     average of the standard normal distribution
     * \param   sigma       sigma value of the standard normal distribution
     * \param   xmin        minimal x value
     * \param   xmax        maximal x value
     * \param   cutoff      use a cutoff
     * \return  standard normal random number
     */

    double sndrn(double rnd, double average, double sigma, double xmin, double xmax, bool cutoff);

//----------------------------------------------------------------------------//
    /**
     * \brief function describes standard normal distribution - interface to Integral
     *
     * \param   x   x-value
     * \return  normalized standard normal distribution with sigma = 1
     */

    double function(double x);

//----------------------------------------------------------------------------//

    /**
     * 1d parametrization - interface to Interpolate
     *
     * \param   x   x-value
     * \return  return the standard normal distribution value of x
     */

    double functionInt(double x);

//----------------------------------------------------------------------------//
    // Getter

    bool get_jt()
    {
        return jt;
    }

    Interpolate get_interpolateJ();

//----------------------------------------------------------------------------//
    // Setter

    void set_jt(bool newJT);

};


#endif /* STANDARDNORMAL_H_ */
