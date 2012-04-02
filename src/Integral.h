/*
 * Integral.h
 *
 *  Created on: 02.08.2010
 *      Author: koehne
 */

#ifndef INTEGRAL_H_
#define INTEGRAL_H_

#include "FunctionOfx.h"
#include <vector>
#include "MathModel.h"

// #include "Output.h"

/**
 * This class provides routines for function integration using Romberg method.
 * Include the function to be integrated in a class that implements the interface FunctionOfx (defined below).
 * Methods contained here are based on the Numerical Recipes (W. H. Press et al.) with the following modifications:
 * <ul><li>Algorithm of route choice in the interpolate method is different: it keeps the partial approximations
 * centered on the target x by keeping two of the grid points closest to x in the middle as long as possible.</li>
 * <li>Power of Substitution (POS) is any real number. If it is zero, no substitution is made. Otherwise the
 * substitution is 1/x^POS if POS&gt;0 or 1/(-x)^(-POS), if POS&lt;0. If one of the limits of integration is zero or
 * has opposite sign than that of POS, it is replaced with infinity. If both limits of integration are zero or have
 * opposite sign than that of POS, integral is considered to have equal to each other (and to infinity) limits, and
 * the returned value is zero.</li></ul>
 * It is possible to evaluate x(rand) such that the integral from xmin to x(rand) is a fraction rand of the original
 * full integral from xmin to xmax. Set the ratio 0&lt;rand&lt;1 as the last argument of one of the open integration
 * methods (integrateOpened or integrateWithSubstitution), and get x(rand) by the subsequent call to getUpperLimit.
 * If rand&lt;0, then -rand is assumed to be the absolute value of the portion of the original integral such that the
 * integral from xmin to x(rand) is equal to this portion. Its sign is determined as the sign of the whole integral
 * (or, rather, the N-1st approximation to its value). If rand is given, it is generally assumed that the function
 * does not change sign on the integration interval. Otherwise, the resulting x(rand) is less predictable. In any
 * case, an approximation to x(rand) is found during evaluation of the original integral, and then refined by the
 * combination of the Newton-Raphson method and bisection.
 * <pre>
 * interface FunctionOfx{
 *     double function(double x);
 * }
 * </pre>
 * @author Dmitry Chirkin
 */

class Integral: public MathModel{

private:

	double max, min;
	double precision;
	int maxSteps;
	int romberg;

	double powerOfSubstitution;
	FunctionOfx *function2use;

	double integralValue;
	double integralError;
	std::vector<double> iX;
	std::vector<double> iY;

	std::vector<double> c;
	std::vector<double> d;
	int romberg4refine;		// set to 2 in constructor

	double randomNumber;
	double randomX;

	bool reverse;			// set to false in Constructor
	double reverseX;
	double savedResult;
	bool randomDo;			// set to false in Constructor
	bool useLog;			// set to false in Constructor

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes class with default settings
     */

public:
	Integral();

       // virtual ~Integral();

    //----------------------------------------------------------------------------------------------------//

    /*!
      initializes class - this is the main constructor, definition of the integration limits, setting of constants
      \f$romberg4refine\f$, \f$pOS\f$, definition whether log-integration or integration with substitution is used
     */

	Integral(int romberg, int maxSteps, double precision);

    //----------------------------------------------------------------------------------------------------//

    /*!
      table of substitutions
     */

	double function(double x);

    //----------------------------------------------------------------------------------------------------//

    /*!
      returns corrected integral by composite trapezoid rule for closed intervals, n=1, 2, 4, 8, ...
      \f[return= \frac{1}{2}[oldSum+resultSum \cdot stepSize]\f]
      formula depends on:
      \f[ \int_a^b f(x)dx=\frac{b-a}{n} [\frac{f(a)+f(b)}{2}+ \sum_{k=1}^{n-1}f(a+k\frac{b-a}{n})] \f]
     */

	double trapezoid(int n, double oldSum);

    //----------------------------------------------------------------------------------------------------//

    /*!
      returns corrected integral by composite trapezoid rule for opened intervals, n=1, 3, 9, ...
      \f[ return= \frac{1}{3} oldSum +resultSum \cdot stepSize; \f]
      formula depends on:
      \f[ \int_a^b f(x)dx=\frac{b-a}{n} \sum_{i=1}^{n}f(a+k\frac{b-a}{n}(\frac{2i-1}{2})) \f]
     */

	double trapezoid3(int n, double oldSum);

    //----------------------------------------------------------------------------------------------------//

    /*!
      returns corrected integral by opened trapezoid rule, n=1, 3, 9, ...
      and computes the approximation to the value of the x(rand)
     */

	double trapezoid3S(int n, double oldSum, int stepNumber);

    //----------------------------------------------------------------------------------------------------//

    /*!
      finds f(x) for f: iY[i]=f(iX[i]), i=start, ..., start+romberg
     */

	void interpolate(int start, double x);

    //----------------------------------------------------------------------------------------------------//

    /*!
      finds integral for closed intervals calculates integral value using trapezoid, k is number of steps;
     */

	double rombergIntegrateClosed();

    //----------------------------------------------------------------------------------------------------//

    /*!
      finds integral for opened intervals, using trapezoid3;
     */

	double rombergIntegrateOpened();

    //----------------------------------------------------------------------------------------------------//

    /*!
      finds integral for opened intervals analog to rombergIntegrateOpened; precision of the result
      is evaluated with respect to the value provided in the argument
     */

	double rombergIntegrateOpened(double bigValue);

    //----------------------------------------------------------------------------------------------------//

    /*!
      finds integral for closed intervals depending on rombergIntegrateClosed;
      \f[ return=aux \cdot rombergClosed \f]
      with \f$ aux=1 \f$ or \f$ aux=-1\f$, depends on the integration limits (\f$min<max\f$ or \f$min>max\f$)
     */

	double integrateClosed(double min, double max, FunctionOfx *function2use);

    //----------------------------------------------------------------------------------------------------//

    /*!
      finds integral for opened intervals;
      analog to integrateClosed, depending on rombergIntegrateOpened
     */

	double integrateOpened(double min, double max, FunctionOfx *function2use);

    //----------------------------------------------------------------------------------------------------//

    /*!
      finds integral for opened intervals
      and computes the value of the x(rand)
     */

	double integrateOpened(double min, double max, FunctionOfx *function2use, double randomRatio);

    //----------------------------------------------------------------------------------------------------//

    /*!
      finds integral for opened intervals using substitution x -&gt; 1/x^(powerOfSubstitution)
     */

	double integrateWithSubstitution(double min, double max, FunctionOfx *function2use, double powerOfSubstitution);

    //----------------------------------------------------------------------------------------------------//

    /*!
      finds integral for opened intervals using substitution x -&gt; 1/x^(powerOfSubstitution)
      and computes the value of the x(rand)
     */

	double integrateWithSubstitution(double min, double max, FunctionOfx *function2use, double powerOfSubstitution, double randomRatio);

    //----------------------------------------------------------------------------------------------------//

    /*!
      using Newton's method refines the value of the upper limit
      that results in the ratio of integrals equal to randomNumber
     */

	void refineUpperLimit(double result);

    //----------------------------------------------------------------------------------------------------//

    /*!
      refines and returns the value of the upper limit x(rand)
     */

	double getUpperLimit();

    //----------------------------------------------------------------------------------------------------//

    /*!
      finds integral for opened intervals using log substitution; \f$ x=\ln(x)\f$
     */

	double integrateWithLog(double min, double max, FunctionOfx *function2use);

    //----------------------------------------------------------------------------------------------------//

    /*!
      finds integral for opened intervals
      and computes the value of the x(rand)
     */

	double integrateWithLog(double min, double max, FunctionOfx *function2use, double randomRatio);

    //----------------------------------------------------------------------------------------------------//

    /*!
      finds integral for opened intervals using substitution x -&gt; 1/log(x)^(powerOfSubstitution)
     */

	double integrateWithLogSubstitution(double min, double max, FunctionOfx *function2use, double powerOfSubstitution);

	// Getters

	double getmaxSteps () { return maxSteps; }

        // Setters

        void set_function2use (FunctionOfx *function) {
            std::cerr << "Integral::set_funtion2use is depricated and might even be buggy. \n better make use of Integral::integrateOpened(...) to set the function to use, and its range...\n";
            this->function2use = function;
        };
};

#endif /* INTEGRAL_H_ */
