/*! \file   StandardNormal.h
*   \brief  Header file for the standard normal routines.
*
*   For more details see the class documentation.
*
*   \date   18.04.2013
*   \author Jan-Hendrik Koehne
*/


#ifndef STANDARDNORMAL_H_
#define STANDARDNORMAL_H_

#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Interpolant.h"
#include <cmath>

/**
 * \brief This class provides routines for evaluation of
 * random numbers distributed normally.
 */

class StandardNormal
{

protected:
    Integral    *integral_;
    Interpolant *interpolant_;

    double      val1_;
    double      val2_;

    bool        do_interpolation_;
    double      norm_;
    int         order_of_interpolation_;
//----------------------------------------------------------------------------//
    //Memberfunctions
    /**
     * \brief class initializer
     *
     * \param   romberg     order of interpolation and integration
     * \param   maxSteps    number of sampling points
     * \param   precision   integration precision
     */

    void Init(int romberg, int maxSteps, double precision);

//----------------------------------------------------------------------------//
    /**
     * \breif evaluates the integrated probability
     *
     * \param   x   wanted value
     * \return  integrated probability
     */

    double IntegratedProbability(double x);
//----------------------------------------------------------------------------//
    /**
     * evaluates the standard normal random number
     *
     * \param   x   wanted value
     * \return  standard normal random number
     */

    double StandardNormalRandomNumber(double x);

//----------------------------------------------------------------------------//

public:

//----------------------------------------------------------------------------//
    //Constructors

    /**
     * \brief initializes class with default settings
     */

    StandardNormal();
    StandardNormal(const StandardNormal&);
    StandardNormal& operator=(const StandardNormal&);

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

//    double sndrn(double rnd, double average, double sigma, double xmin, double xmax, bool cutoff);
    double StandardNormalRandomNumber(double rnd, double average, double sigma, double xmin, double xmax, bool cutoff);

//----------------------------------------------------------------------------//
    /**
     * \brief function describes standard normal distribution - interface to Integral
     *
     * \param   x   x-value
     * \return  normalized standard normal distribution with sigma = 1
     */

    double FunctionToIntegral(double x);

//----------------------------------------------------------------------------//

    /**
     * 1d parametrization - interface to Interpolate
     *
     * \param   x   x-value
     * \return  return the standard normal distribution value of x
     */

    double FunctionToBuildInterpolant(double x);

//----------------------------------------------------------------------------//

    void EnableInterpolation();

//----------------------------------------------------------------------------//

    void DisableInterpolation();

//----------------------------------------------------------------------------//
    // Getter


//----------------------------------------------------------------------------//
    // Setter


};


#endif /* STANDARDNORMAL_H_ */
