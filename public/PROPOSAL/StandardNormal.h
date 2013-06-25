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
    Integral*       integral_;
    Interpolant*    interpolant_;

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

public:

//----------------------------------------------------------------------------//
    //Constructors

    /**
     * \brief initializes class with default settings
     */

    StandardNormal();
    StandardNormal(const StandardNormal&);
    StandardNormal& operator=(const StandardNormal&);
    bool operator==(const StandardNormal &normal) const;
    bool operator!=(const StandardNormal &normal) const;

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

    void swap(StandardNormal &normal);
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

    double StandardNormalRandomNumber(double rnd, double average, double sigma, double xmin, double xmax, bool cutoff);

//----------------------------------------------------------------------------//

    void EnableInterpolation(std::string path ="", bool raw=false);

//----------------------------------------------------------------------------//

    void DisableInterpolation();

//----------------------------------------------------------------------------//
    //Getter
	bool GetDoInterpolation() const {
		return do_interpolation_;
	}

	Integral* GetIntegral() const {
		return integral_;
	}

	Interpolant* GetInterpolant() const {
		return interpolant_;
	}

	double GetNorm() const {
		return norm_;
	}

	int GetOrderOfInterpolation() const {
		return order_of_interpolation_;
	}

	double GetVal1() const {
		return val1_;
	}

	double GetVal2() const {
		return val2_;
	}
//----------------------------------------------------------------------------//
    //Setter
	void SetDoInterpolation(bool doInterpolation);
	void SetIntegral(Integral* integral);
	void SetInterpolant(Interpolant* interpolant);
	void SetNorm(double norm);
	void SetOrderOfInterpolation(int orderOfInterpolation);
	void SetVal1(double val1);
	void SetVal2(double val2);
};


#endif /* STANDARDNORMAL_H_ */
