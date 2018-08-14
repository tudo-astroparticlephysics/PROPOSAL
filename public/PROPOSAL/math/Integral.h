
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/


#pragma once

#include <functional>
#include <iostream>
#include <vector>

namespace PROPOSAL {

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

class Integral
{

public:
    /**
     * initializes class with default settings
     */
    Integral();

    Integral(const Integral&);

    /*!
     * initializes class - this is the main constructor, definition of the
     * integration limits, setting of constants \f$romberg4refine\f$, \f$pOS\f$,
     * definition whether log-integration or integration with
     * substitution is used
     */
    Integral(int romberg, int maxSteps, double precision);

    ~Integral();

    Integral& operator=(const Integral&);
    bool operator==(const Integral& integral) const;
    bool operator!=(const Integral& integral) const;

    void swap(Integral& integral);

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

    double Integrate(double min,
                     double max,
                     std::function<double(double)> integrand,
                     int method,
                     double powerOfSubstitution = 0);

    //----------------------------------------------------------------------------//

    /*!
     * finds integral: choose the in integration method with the last parameter
     * like in Integrate but just for two cases
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
     * \param   randomRatio     Random Ratio in which the old sum is weighted
     * \return  Integration result
     */

    double IntegrateWithRandomRatio(double min,
                                    double max,
                                    std::function<double(double)> integrand,
                                    int method,
                                    double randomRatio,
                                    double powerOfSubstitution = 0);

    //----------------------------------------------------------------------------//

    // --------------------------------------------------------------------- //
    // Getter
    // --------------------------------------------------------------------- //

    std::function<double(double)> GetIntegrand() const { return integrand_; }
    double GetMax() const { return max_; }
    double GetMin() const { return min_; }
    double GetPowerOfSubstitution() const { return powerOfSubstitution_; }
    double GetPrecision() const { return precision_; }
    double GetRandomNumber() const { return randomNumber_; }
    double GetRandomX() const { return randomX_; }
    double GetReverseX() const { return reverseX_; }
    double GetSavedResult() const { return savedResult_; }
    int GetMaxStepsUpperLimit() const { return maxSteps_upper_limit_; }
    int GetMaxStepsRomberg() const { return maxSteps_romberg_; }
    int GetRomberg() const { return romberg_; }
    int GetRomberg4refine() const { return romberg4refine_; }
    bool GetRandomDo() const { return randomDo_; }
    bool GetReverse() const { return reverse_; }
    bool GetUseLog() const { return useLog_; }

    // --------------------------------------------------------------------- //
    // Setter
    // --------------------------------------------------------------------- //

    void Set_Function(std::function<double(double)> integrand)
    {
        std::cerr << "Integral::set_funtion2use is depricated and might even be buggy. \n better make use of "
                     "Integral::integrateOpened(...) to set the function to use, and its range...\n";
        this->integrand_ = integrand;
    }
    void SetIntegrand(std::function<double(double)> integrand);
    void SetMax(double max);
    void SetMaxStepsUpperLimit(int maxSteps);
    void SetMaxStepsRomberg(int maxSteps);
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

private:
    struct InterpolationResults
    {
        double Value;
        double Error;
    };

    struct QuadpackResults
    {
        QuadpackResults()
            : value(0.0)
            , abserr(0.0)
            , neval(0)
            , ier(0)
        {
        }

        double value;
        double abserr;
        int neval;
        int ier;
    };

private:
    int maxSteps_romberg_;
    int maxSteps_upper_limit_;
    int romberg_;
    double precision_;
    double max_, min_;

    std::vector<double> iX_;
    std::vector<double> iY_;

    std::vector<double> c_;
    std::vector<double> d_;

    std::function<double(double)> integrand_;

    int romberg4refine_; // set to 2 in constructor
    double powerOfSubstitution_;
    bool randomDo_; // set to false in Constructor
    bool useLog_;   // set to false in Constructor
    double randomNumber_;
    double randomX_;
    bool reverse_; // set to false in Constructor
    double reverseX_;
    double savedResult_;

    std::vector<double> q_last_3_results_;
    std::vector<double> q_rlist2_; // epstab
    std::vector<double> q_iord_;

    // ----------------------------------------------------------------------------
    /// @brief This function is a translation of the fortran 77 subroutine
    ///        dqk21 from the package QUADPACK by Piessens et al. (1983),
    ///        which is needed by qags.
    // ----------------------------------------------------------------------------
    std::pair<Integral::InterpolationResults, Integral::InterpolationResults> q_gaus_kronrod_21(double q_min,
                                                                                                double q_max);

    // ----------------------------------------------------------------------------
    /// @brief This function is a translation of the fortran 77 subroutine
    ///        dqelg from the package QUADPACK by Piessens et al. (1983),
    ///        which is needed by qags.
    // ----------------------------------------------------------------------------
    Integral::InterpolationResults q_epsilon_extrapolation(int q_limit_epsilon_table, int numrl2, int nres);

    // ----------------------------------------------------------------------------
    /// @brief This function is a translation of the fortran 77 subroutine
    ///        dqpsrt from the package QUADPACK by Piessens et al. (1983),
    ///        which is needed by qags.
    // ----------------------------------------------------------------------------
    int q_sort(int q_limit, int last, int maxerr, int nrmax, const std::vector<double>& q_elist);

    // ----------------------------------------------------------------------------
    /// @brief QUADPACK implementation of the gauss kronrod integration.
    ///        This function is a translation of the fortran 77 subroutine
    ///        qags from the package QUADPACK by Piessens et al. (1983).
    ///
    /// @param limit        determines the maximum number of subintervals
    ///                     in the partition of the given integration interval
    ///                     (a,b), limit > 1.
    /// @param q_epsabs_    absolute accoracy requested
    /// @param q_epsrel_    relative accuracy requested
    ///
    /// @return             QuadpackResults, a struct containing
    //
    ///                     value: The approximation of the integral
    ///                     abserr: The approximation of the absolute error
    ///                     neval: Number of function evaluations
    ///                     ier: The error code
    ///
    ///                     ier = 0 normal and reliable termination of the
    ///                             routine. it is assumed that the requested
    ///                             accuracy has been achieved.
    ///                     ier > 0 abnormal termination of the routine
    ///                             the estimates for result and error are
    ///                             less reliable. it is assumed that the
    ///                             requested accuracy has not been achieved.
    ///
    ///                     error messages
    ///                     ier = 1 maximum number of subdivisions allowed
    ///                             has been achieved. one can allow more
    ///                             subdivisions by increasing the value of
    ///                             limit (and taking the according dimension
    ///                             adjustments into account). however, if
    ///                             this yield no improvement it is advised
    ///                             to analyze the integrand in order to
    ///                             determine the integration difficulaties.
    ///                             if the position of a local difficulty can
    ///                             be determined (i.e.singularity,
    ///                             discontinuity within the interval) one
    ///                             will probably gain from splitting up the
    ///                             interval at this point and calling the
    ///                             integrator on the subranges. if possible,
    ///                             an appropriate special-purpose integrator
    ///                             should be used which is designed for
    ///                             handling the type of difficulty involved.
    ///                         = 2 the occurrence of roundoff error is
    ///                             detected, which prevents the requested
    ///                             tolerance from being achieved.
    ///                         = 3 extremely bad integrand behaviour occurs
    ///                             at some points of the integration
    ///                             interval.
    ///                         = 6 the input is invalid, because
    ///                             (epsabs <= 0 and epsrel <= 0
    ///                             or limit < 1).
    ///                             result, abserr, neval, last are set
    ///                             to zero.
    // ----------------------------------------------------------------------------
    QuadpackResults qags(double limit = 50, double q_epsabs = 1.0e-50, double q_epsrel = 1.0e-6);

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

    double InitIntegralOpenedAndClosed(double min, double max, std::function<double(double)> integrand);

    //----------------------------------------------------------------------------//

    double InitIntegralWithSubstitution(double min,
                                        double max,
                                        std::function<double(double)> integrand,
                                        double powerOfSubstitution);

    //----------------------------------------------------------------------------//

    double InitIntegralWithLogSubstitution(double min,
                                           double max,
                                           std::function<double(double)> integrand,
                                           double powerOfSubstitution);

    //----------------------------------------------------------------------------//

    double InitIntegralWithLog(double min, double max, std::function<double(double)> integrand);

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
    double IntegrateClosed(double min, double max, std::function<double(double)> integrand);

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
    double IntegrateOpened(double min, double max, std::function<double(double)> integrand);

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

    double IntegrateWithSubstitution(double min,
                                     double max,
                                     std::function<double(double)> integrand,
                                     double powerOfSubstitution);

    //----------------------------------------------------------------------------//

    /*!
     * finds integral for opened intervals using log substitution; \f$ x=\ln(x)\f$
     *
     * \param   min             lower integration limit
     * \param   max             upper integration limit
     * \param   function2use    integrand
     * \return  Integration result
     */

    double IntegrateWithLog(double min, double max, std::function<double(double)> integrand);

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

    double IntegrateWithLogSubstitution(double min,
                                        double max,
                                        std::function<double(double)> integrand,
                                        double powerOfSubstitution);

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

    double IntegrateWithLog(double min, double max, std::function<double(double)> integrand, double randomRatio);

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

    double IntegrateWithSubstitution(double min,
                                     double max,
                                     std::function<double(double)> integrand,
                                     double powerOfSubstitution,
                                     double randomRatio);
};

} // namespace PROPOSAL
