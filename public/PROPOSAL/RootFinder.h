/**
 * Original comment of the Java class by Dmitry Chirkin:
 * This class provides routines for calculating roots by
 * the combination of the Newton-Raphson method and bisection.
 * Include the function to be integrated in a class that implements
 * the interface FunctionOfx (defined below).
 * Methods contained here are based on the Numerical Recipes
 * (W. H. Press et al.).
 *
 * @author Dmitry Chirkin
 * Ported to C++ Jan-Hendrik Köhne
 */

#ifndef FINDROOT_H
#define FINDROOT_H

#include "boost/function.hpp"

class FindRoot
{
private:

    double  _precision;
    int     _maxSteps;

    boost::function<double (double)> function_;
    boost::function<double (double)> differentiated_function_;

//----------------------------------------------------------------------------//
    //Memberfunctions

    double function(double x);

//----------------------------------------------------------------------------//

    double dFunction(double x);

//----------------------------------------------------------------------------//

public:

    //Constructors

    FindRoot();
    FindRoot(const FindRoot&);
    FindRoot& operator=(const FindRoot&);
	FindRoot(int maxSteps, double precision);

//----------------------------------------------------------------------------//
    //Destructor

    ~FindRoot();
//----------------------------------------------------------------------------//
    //Memberfunctions

	double findRoot(double min, double max, double startX, DFunctionOfx *function2use, double rightSide);


  


};

#endif // _FINDROOT_H
