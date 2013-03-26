/*! \file   Interpolant.h
*   \brief  Header file for definition of the Interpolant class object.
*
*   For more details see the class documentation.
*
*   \date   02.08.2010
*   \author Jan-Hendrik Koehne
*/

#ifndef INTERPOLANT_H_
#define INTERPOLANT_H_

#include <vector>
#include <cmath>

#include "PROPOSAL/MathModel.h"
#include "PROPOSAL/FunctionInt2.h"
#include "PROPOSAL/FunctionInt.h"
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include "Output.h"


/**
 *\class Interpolant
 *
 * This class provides routines for function interpolation. Include the function to be interpolated
 * in a class that implements the interface FunctionInt or FunctionInt2 (defined below).
 * Methods contained here are based on the Numerical Recipes (W. H. Press et al.).
 * <pre>
 * interface FunctionInt{
 *     double functionInt(double x);
 * }
 *
 * interface FunctionInt2{
 *     double functionInt(double x1, double x2);
 * }
 * </pre>
 * @author Dmitry Chirkin
 */

class Interpolant{ /// implements FunctionInt{

private:
    int     romberg_, rombergY_;
    double  *iX_;
    double  *iY_;
    double  *c_;
    double  *d_;
    int     max_;
    double  xmin_, xmax_, step_;
    bool    rational_, relative_;

    //FunctionInt2*   function2int;
    boost::function<double (double)> function1d_;
    boost::function<double (double,double)> function2d_;
    Interpolant*    Interpolant_;

    int     row_, starti_;
    bool    rationalY_, relativeY_;

    bool    reverse_, self_, flag_;		// Self is setted to true in constructor
    bool    isLog_, logSubst_;

    double  precision_, worstX_;
    double  precision2_, worstX2_;
    double  precisionY_, worstY_;

    bool    fast_;			// Is setted to true in constructor

    double  x_save_, y_save_;	// Is setted to 1 and 0 in constructor

    Output  *output_;

    /*!
    * interpolates f(x) based on the values iY[i]=f(iX[i]) in the romberg-vicinity of x
    *
    * \param   x        position of the function
    * \param   start    start position of the sampling points for interpolation
    * \return  Interpolation result
    */
    double interpolate(double x, int start);

    //------------------------------------------------------------------------//



protected:

    const static double bigNumber_;
    const static double aBigNumber_;


    //------------------------------------------------------------------------//


public:

    /*!
    * Default Constructor.
    * Default Constructor, which does nothing.
    */
    Interpolant();

    /*!
    * Main constructor for the 1-dimensional functions.
    *
    * Main constructor and initializes class for the 1-dimensional function.
    *
    * \param   max          number of sampling points
    * \param   xmin         lower limit of the interpolation table
    * \param   xmax         upper limit of the interpolation table
    * \param   function2int function which will be interpolated
    * \param   romberg      Order of interpolation
    * \param   rational     interpolate with rational function?
    * \param   relative     save error relative to the function value?
    * \param   isLog        substitute x = log(x)?
    * \param   rombergY     Order of interpolation for inverse interpolation(Find Limit)
    * \param   rationalY    interpolate with rational function?
    * \param   relativeY    save error relative to the interpolated x-value?
    * \param   logSubst     substitute f(x) = log(f(x))?
    * \return
    */
    Interpolant(int max, double xmin, double xmax, boost::function<double (double)> function1d,
                int romberg, bool rational, bool relative, bool isLog,
                int rombergY, bool rationalY, bool relativeY, bool logSubst);

    //------------------------------------------------------------------------//

    /*!
    * Main constructor for the 2-dimensional functions.
    *
    * Main constructor and initializes class for the 2-dimensional function.
    *
    * \param   max          number of sampling points
    * \param   xmin         lower limit of the interpolation table
    * \param   xmax         upper limit of the interpolation table
    * \param   function2int function which will be interpolated
    * \param   romberg1     Order of interpolation f(x,a)
    * \param   rational1    interpolate with rational function?
    * \param   relative1    save error relative to the function value?
    * \param   isLog1       substitute x1 = log(x2)
    * \param   romberg2     Order of interpolation for f(a,x)
    * \param   rational2    interpolate with rational function?
    * \param   relative2    save error relative to the function value?
    * \param   isLog2       substitute x2 = log(x2)
    * \param   rombergY     Order of interpolation for inverse interpolation(Find Limit)
    * \param   rationalY    interpolate with rational function?
    * \param   relativeY    save error relative to the interpolated x-value
    * \param   logSubst     substitute f(x1,x2) = log(f(x1,x2))
    * \return
    */
    Interpolant(int max1, double x1min, double x1max, int max2, double x2min, double x2max,  boost::function<double (double,double)> function2d,
                int romberg1, bool rational1, bool relative1, bool isLog1,
                int romberg2, bool rational2, bool relative2, bool isLog2,
                int rombergY, bool rationalY, bool relativeY, bool logSubst);

    //------------------------------------------------------------------------//


    /*!
    * Constructor for the 1-dimensional functions if the array already exists.
    *
    * Constructs and initializes class for the 1-dimensional function with the
    * given array of x and f(x) values.
    *
    * \param   x    vector of x values.
    * \param   xmin vector of f(x) values
    * \param   romberg      Order of interpolation
    * \param   rational     interpolate with rational function
    * \param   relative     save error relative to the function value?
    * \return
    */
    Interpolant(std::vector<double> x, std::vector<double> y, int romberg, bool rational, bool relative);

    //------------------------------------------------------------------------//

    /**
    * Auxiliary class initializer.
    *
    * Initializes class for the 1-dimensional function.
    *
    * \param   max          number of sampling points
    * \param   xmin         lower limit of the interpolation table
    * \param   xmax         upper limit of the interpolation table
    * \param   romberg      order of interpolation
    * \param   rational     interpolate with rational function
    * \param   relative     save error relative to the function value
    * \param   isLog        substitute x = log(x)
    * \param   rombergY     order of interpolation for inverse interpolation(Find Limit)
    * \param   rationalY    interpolate with rational function
    * \param   relativeY    save error relative to the interpolated x-value
    * \param   logSubst     substitute f(x) = log(f(x))
    * \return
    */
    void InitInterpolant(int max, double xmin, double xmax,
                         int romberg, bool rational, bool relative, bool isLog,
                         int rombergY, bool rationalY, bool relativeY, bool logSubst);

    //------------------------------------------------------------------------//

    /**
    * Defines a function for every row.
    *
    * \param    x
    * \return   Function value f(x)
    */

    double functionInt(double x);

    //------------------------------------------------------------------------//

    /**
    * Interpolates f(x) for 1d function
    *
    * \param    x
    * \return   interpolated value f(x)
    */

    double interpolate(double x);

    //------------------------------------------------------------------------//

    /**
    * Interpolates f(x) for 2d function
    *
    * \param    x1
    * \param    x2
    * \return   interpolated value f(x1,x2)
    */

    double interpolate(double x1, double x2);

    //------------------------------------------------------------------------//

    /**
    * Interpolates f(x) for 1d function if the arrays already exist.
    *
    * \param    x
    * \return   interpolated value f(x)
    */

    double InterpolateArray(double x);

    //------------------------------------------------------------------------//



    /**
    * Finds x: f(x)=y., 1d initialization required
    *
    * Finds x: f(x)=y. 1d initialization required.
    *
    * \param    y
    * \return   interpolated value x(y);
    */

    double findLimit(double y);

    //------------------------------------------------------------------------//

    /**
    * Finds x: f(a,x)=y.
    *
    * Finds x: f(a,x)=y, 2d initialization required
    *
    * \param    x1
    * \param    y
    * \return   interpolated value x(y);
    */

    double findLimit(double x1, double y);

    //------------------------------------------------------------------------//

    double Get2dFunctionFixedY(double x); //TOMSASZ
    //------------------------------------------------------------------------//

    /**
    * Exp(x) with CutOff.
    *
    * if x > exp(aBigNumber): exp(x) \n
    * else: 0
    *
    * \param    x
    * \return   exp(x) OR 0;
    */

    double exp(double x);

    //------------------------------------------------------------------------//

    /**
    * Log if not zero.
    *
    * if x > 0: log(x) \n
    * else: bigNumber
    *
    * \param    x
    * \return   log(x) OR bigNumber;
    */

    double log(double x);

    //------------------------------------------------------------------------//

    /**
    * Log it again.
    *
    * \param    x
    * \return   log(x) OR y_save;
    */

    double slog(double x);

    //------------------------------------------------------------------------//

    /**
    * Getter for Interpolant object.
    *
    * \return   Interpolant object;
    */
    Interpolant* get_Interpolant()
    {
        return Interpolant_;
    }

    /**
    * Getter for max.
    *
    * \return   the number of sampling points;
    */
    int get_max()
    {
        return max_;
    }
};

#endif /* INTERPOLANT_H_ */
