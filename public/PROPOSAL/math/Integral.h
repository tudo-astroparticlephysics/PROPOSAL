/*! \file   Integral.h
*   \brief  Header file for definition of the integral class object.
*
*   For more details see the class documentation.
*
*   \date   02.08.2010
*   \author Jan-Hendrik Koehne
*/
#pragma once

#ifndef INTEGRAL_H_
#define INTEGRAL_H_

#include <vector>
#include <iostream>
#include <boost/function.hpp>


namespace PROPOSAL{

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
 *
 * @author Dmitry Chirkin
 */

class Integral{

private:
    struct InterpolationResults
    {
        double Value;
        double Error;
    };

private:

    int     maxSteps_;
    int     romberg_;
    double  precision_;
    double  max_, min_;


    // double  integralValue_;
    // double  integralError_;

    std::vector<double> iX_;
    std::vector<double> iY_;

    std::vector<double> c_;
    std::vector<double> d_;

    boost::function<double (double)> integrand_;


    int     romberg4refine_;		// set to 2 in constructor
    double  powerOfSubstitution_;
    bool    randomDo_;			// set to false in Constructor
    bool    useLog_;			// set to false in Constructor
    double  randomNumber_;
    double  randomX_;
    bool    reverse_;			// set to false in Constructor
    double  reverseX_;
    double  savedResult_;

    int q_limit_;
    // double q_epsabs_;
    // double q_epsrel_;
    std::vector<double> q_last_3_results_;
    std::vector<double> q_rlist2_; //epstab
    std::vector<double> q_iord_;
    std::vector<double> q_alist_;
    std::vector<double> q_blist_;
    std::vector<double> q_elist_;
    std::vector<double> q_rlist_;
    std::vector<double> q_fv1_;
    std::vector<double> q_fv2_;
    // const double q_weights_gaus_10p_ [5] =
    // {
    //     0.066671344308688137593568809893332,
    //     0.149451349150580593145776339657697,
    //     0.219086362515982043995534934228163,
    //     0.269266719309996355091226921569469,
    //     0.295524224714752870173892994651338
    // };
    // const double q_weights_kronrod_21p_ [11] =
    // {
    //     0.011694638867371874278064396062192,
    //     0.032558162307964727478818972459390,
    //     0.054755896574351996031381300244580,
    //     0.075039674810919952767043140916190,
    //     0.093125454583697605535065465083366,
    //     0.109387158802297641899210590325805,
    //     0.123491976262065851077958109831074,
    //     0.134709217311473325928054001771707,
    //     0.142775938577060080797094273138717,
    //     0.147739104901338491374841515972068,
    //     0.149445554002916905664936468389821
    // };
    // const double q_abscissae_kronrod_21p_[11] =
    // {
    //     0.995657163025808080735527280689003,
    //     0.973906528517171720077964012084452,
    //     0.930157491355708226001207180059508,
    //     0.865063366688984510732096688423493,
    //     0.780817726586416897063717578345042,
    //     0.679409568299024406234327365114874,
    //     0.562757134668604683339000099272694,
    //     0.433395394129247190799265943165784,
    //     0.294392862701460198131126603103866,
    //     0.148874338981631210884826001129720,
    //     0.000000000000000000000000000000000
    // };
    // const double q_epmach_ = std::numeric_limits<double>::epsilon(); // machine epsilon
    // const double q_uflow_ = std::numeric_limits<double>::min(); // smallest finite value
    // const double q_oflow_ = std::numeric_limits<double>::max(); // largest finite value
    // limexp is the maximum number of elements the epsilon table can contain.
    // if this number is reached, the upper diagonal of the epsilon table is deleted.
    // const int q_limit_epsilon_table_ = 50;

    std::pair<Integral::InterpolationResults, Integral::InterpolationResults> q_gaus_kronrod_21(double q_min, double q_max);
    Integral::InterpolationResults q_epsilon_extrapolation(int numrl2, int nres);
    int q_sort(int last, int maxerr, int nrmax);
    double qags();


//----------------------------------------------------------------------------//

    /*!
     * returns corrected integral by composite trapezoid rule for
     * closed intervals, n=1, 2, 4, 8, ...
     * \f[return= \frac{1}{2}[oldSum+resultSum \cdot stepSize]\f]
     * formula depends on:
     * \f[ \int_a^b f(x)dx=\frac{b-a}{n} [\frac{f(a)+f(b)}{2}
     * + \sum_{k=1}^{n-1}f(a+k\frac{b- a}{n})] \f]
     *
     * \param   n         Number of sampling points
     * \param   oldSum    Old integration result
     * \return  Integration Result
     */
    double Trapezoid(int n, double oldSum);

//----------------------------------------------------------------------------//

    /*!
     * returns corrected integral by composite trapezoid
     * rule for opened intervals, n=1, 3, 9, ...
     * \f[ return= \frac{1}{3} oldSum +resultSum \cdot stepSize; \f]
     * formula depends on:
     * \f[ \int_a^b f(x)dx=\frac{b-a}{n} \sum_{i=1}^{n}f(a+k\frac{b-a}{n}(\frac{2i-1}{2})) \f]
     *
     * \param   n           Number of sampling points
     * \param   oldSum      Old integration result
     * \return  Integration Result
     */
   double Trapezoid3(int n, double oldSum);

//----------------------------------------------------------------------------//

   /*!
    * returns corrected integral by opened trapezoid rule, n=1, 3, 9, ...
    * and computes the approximation to the value of the x(rand)
    *
    * \param   n           Number of sampling points
    * \param   oldSum      Old integration result
    * \return  Integration result
    */
   double Trapezoid3S(int n, double oldSum, int stepNumber);

//----------------------------------------------------------------------------//

   /*!
    * Interface for the integrand. This function does all the substitutions (like
    * log, exp, pow...), which are needed by the integration routines
    * (e.g. Trapezoid).
    * \param   x
    * \return  modified integrand value
    */
   double Function(double x);

//----------------------------------------------------------------------------//


   /*!
    * Interpolates the integral value by interpolating/extrapolating the value
    * on the basis of the last "romberg" integral values. Mostly it is used with
    * x=0 because it extrapolates to an infinite number of sampling points or
    * in other words a step size of 0. So it extrapolates the function: integral
    * values over stepsize(or 1/sampling points).
    *
    * \param   start   Starting point of the interpolation
    * \param   x       Searched function value for f(x)
    * \return  Interpolated result and error as struct for f(x)
    */
   Integral::InterpolationResults Interpolate(int start, double x);

//----------------------------------------------------------------------------//


   /*!
    * finds integral for closed intervals calculates integral value
    * using trapezoid, k is number of steps;
    *
    * \return  Integration result
    */
   double RombergIntegrateClosed();

//----------------------------------------------------------------------------//

   /*!
    * finds integral for opened intervals, using trapezoid3;
    *
    * \return  Integration result
    */
   double RombergIntegrateOpened();

//----------------------------------------------------------------------------//

   /*!
    * finds integral for opened intervals analog to rombergIntegrateOpened;
    * precision of the result is evaluated with respect
    * to the value provided in the argument
    *
    * \param   bigValue    integration error has to be < bigValue*precision
    * \return  Integration result
    */
   double RombergIntegrateOpened(double bigValue);

//----------------------------------------------------------------------------//

   /*!
    * using Newton's method refines the value of the upper limit
    * that results in the ratio of integrals equal to randomNumber:
    *
    * \f$ \int_{low}^y f(x) dx = \int_{low}^{hi} f(x) dx * rnd \f$
    *
    * \param   result  \f$\int_{low}^{hi} f(x) dx\f$
    * \return	y
    */

   void RefineUpperLimit(double result);

//----------------------------------------------------------------------------//

   double InitIntegralOpenedAndClosed(double min, double max, boost::function<double (double)> integrand);

//----------------------------------------------------------------------------//

   double InitIntegralWithSubstitution(double min, double max, boost::function<double (double)> integrand, double powerOfSubstitution);

//----------------------------------------------------------------------------//
   double InitIntegralWithLogSubstitution(double min, double max, boost::function<double (double)> integrand, double powerOfSubstitution);

//----------------------------------------------------------------------------//

   double InitIntegralWithLog(double min, double max, boost::function<double (double)> integrand);


//----------------------------------------------------------------------------//

   /*!
    * finds integral for closed intervals depending on rombergIntegrateClosed;
    * \f[ return=aux \cdot rombergClosed \f]
    * with \f$ aux=1 \f$ or \f$ aux=-1\f$, depends on
    * the integration limits (\f$min<max\f$ or \f$min>max\f$)
    *
    * \param   min             lower integration limit
    * \param   max             upper integration limit
    * \param   function2use    integrand
    * \return  Integration result
    */
    double IntegrateClosed(double min, double max, boost::function<double (double)> integrand);

//----------------------------------------------------------------------------//

   /*!
    * finds integral for opened intervals;
    * analog to integrateClosed, depending on rombergIntegrateOpened
    *
    * \param   min             lower integration limit
    * \param   max             upper integration limit
    * \param   function2use    integrand
    * \return  Integration result
    */
   double IntegrateOpened(double min, double max, boost::function<double (double)> integrand);

//----------------------------------------------------------------------------//

   /*!
    * finds integral for opened intervals using substitution x -&gt; 1/x^(powerOfSubstitution)
    *
    * \param   min                 lower integration limit
    * \param   max                 upper integration limit
    * \param   function2use        integrand
    * \param   powerOfSubstitution x' = 1/x^(powerOfSubstitution)
    * \return  Integration result
    */

   double IntegrateWithSubstitution(double min, double max, boost::function<double (double)> integrand, double powerOfSubstitution);

//----------------------------------------------------------------------------//

   /*!
    * finds integral for opened intervals using log substitution; \f$ x=\ln(x)\f$
    *
    * \param   min             lower integration limit
    * \param   max             upper integration limit
    * \param   function2use    integrand
    * \return  Integration result
    */

   double IntegrateWithLog(double min, double max, boost::function<double (double)> integrand);

//----------------------------------------------------------------------------//


   /*!
    * finds integral for opened intervals using substitution x -&gt; 1/log(x)^(powerOfSubstitution)
    *
    * \param   min             lower integration limit
    * \param   max             upper integration limit
    * \param   function2use    integrand
    * \param   powerOfSubstitution x' = 1/log(x)^(powerOfSubstitution)
    * \return  Integration result
    */

   double IntegrateWithLogSubstitution(double min, double max, boost::function<double (double)> integrand, double powerOfSubstitution);

//----------------------------------------------------------------------------//
    /*!
     * finds integral for opened intervals
     * and computes the value of the x(rand)
     *
     * \param   min             lower integration limit
     * \param   max             upper integration limit
     * \param   function2use    integrand
     * \param   randomRatio     Random Ratio in which the old sum is weighted
     * \return  Integration result
     */

    double IntegrateWithLog(double min, double max, boost::function<double (double)> integrand, double randomRatio);

//----------------------------------------------------------------------------//

    /*!
     * finds integral for opened intervals using substitution x -&gt; 1/x^(powerOfSubstitution)
     * and computes the value of the x(rand)
     *
     * \param   min             lower integration limit
     * \param   max             upper integration limit
     * \param   function2use    integrand
     * \param   powerOfSubstitution x' = 1/x^(powerOfSubstitution)
     * \param   randomRatio     Random Ratio in which the old sum is weighted
     * \return  Integration result
     */

    double IntegrateWithSubstitution(double min, double max, boost::function<double (double)> integrand, double powerOfSubstitution, double randomRatio);

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

public:

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

    /**
     * initializes class with default settings
     */
    Integral();

//----------------------------------------------------------------------------//

    Integral(const Integral&);
    Integral& operator=(const Integral&);
    bool operator==(const Integral &integral) const;
    bool operator!=(const Integral &integral) const;

//----------------------------------------------------------------------------//

    /*!
     * initializes class - this is the main constructor, definition of the
     * integration limits, setting of constants \f$romberg4refine\f$, \f$pOS\f$,
     * definition whether log-integration or integration with
     * substitution is used
     */
    Integral(int romberg, int maxSteps, double precision);


//----------------------------------------------------------------------------//

    /*!
     * refines and returns the value of the upper limit x(rand)
     *
     * \return UpperLimit of Integral
     */

    double GetUpperLimit();

//----------------------------------------------------------------------------//

    /*!
     * finds integral: choose the in integration method with the last parameter
     *   method = 1: IntegrateClosed
     *   method = 2: IntegrateOpened
     *   method = 3: IntegrateWithSubstitution
     *   method = 4: IntegrateWithLog
     *   method = 5: IntegrateWithLogSubstitution
     *
     * \param   min             lower integration limit
     * \param   max             upper integration limit
     * \param   function2use    integrand
     * \param   method          integration method
     * \return  Integration result
     */

    double Integrate(double min, double max, boost::function<double (double)> integrand, int method, double powerOfSubstitution=0);

 //----------------------------------------------------------------------------//

    /*!
     * finds integral: choose the in integration method with the last parameter
     * like in Integrate but just for two cases
     *   method = 3: IntegrateWithSubstitution
     *   method = 4: IntegrateWithLog
     *
     * \param   min             lower integration limit
     * \param   max             upper integration limit
     * \param   function2use    integrand
     * \param   method          integration method
     * \param   randomRatio     Random Ratio in which the old sum is weighted
     * \return  Integration result
     */

    double IntegrateWithRandomRatio(double min, double max, boost::function<double (double)> integrand, int method, double randomRatio, double powerOfSubstitution=0);

 //----------------------------------------------------------------------------//

    void swap(Integral &integral);

 //----------------------------------------------------------------------------//

    boost::function<double (double)> GetIntegrand() const
    {
        return integrand_;
    }
//----------------------------------------------------------------------------//
// Setters

    void Set_Function(boost::function<double (double)> integrand)
    {
        std::cerr << "Integral::set_funtion2use is depricated and might even be buggy. \n better make use of Integral::integrateOpened(...) to set the function to use, and its range...\n";
        this->integrand_ = integrand;
    }

/*!
 * Destructor
 */
    ~Integral();

// Getters

	double GetMax() const {
		return max_;
	}

	int GetMaxSteps() const {
		return maxSteps_;
	}

	double GetMin() const {
		return min_;
	}

	double GetPowerOfSubstitution() const {
		return powerOfSubstitution_;
	}

	double GetPrecision() const {
		return precision_;
	}

	bool GetRandomDo() const {
		return randomDo_;
	}

	double GetRandomNumber() const {
		return randomNumber_;
	}

	double GetRandomX() const {
		return randomX_;
	}

	bool GetReverse() const {
		return reverse_;
	}

	double GetReverseX() const {
		return reverseX_;
	}

	int GetRomberg() const {
		return romberg_;
	}

	int GetRomberg4refine() const {
		return romberg4refine_;
	}

	double GetSavedResult() const {
		return savedResult_;
	}

	bool GetUseLog() const {
		return useLog_;
	}

	void SetIntegrand(boost::function<double(double)> integrand);
	void SetMax(double max);
	void SetMaxSteps(int maxSteps);
	void SetMin(double min);
	void SetPowerOfSubstitution(double powerOfSubstitution);
	void SetPrecision(double precision);
	void SetRandomDo(bool randomDo);
	void SetRandomNumber(double randomNumber);
	void SetRandomX(double randomX);
	void SetReverse(bool reverse);
	void SetReverseX(double reverseX);
	void SetRomberg(int romberg);
	void SetRomberg4refine(int romberg4refine);
	void SetSavedResult(double savedResult);
	void SetUseLog(bool useLog);

};

}

#endif /*INTEGRAL_H_ */
