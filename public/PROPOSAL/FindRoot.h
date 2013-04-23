// FindRoot header file


/**
 * Original comment of the Java class by Dmitry Chirkin:
 * This class provides routines for calculating roots by
 * the combination of the Newton-Raphson method and bisection.
 * Include the function to be integrated in a class that implements
 * the interface FunctionOfx (defined below).
 * Methods contained here are based on the Numerical Recipes
 * (W. H. Press et al.).
 * <pre>
 * interface DFunctionOfx extends FunctionOfx{
 *     double dFunction(double x);
 * }
 * </pre>
 * For the definition of interface FunctionOfx
 * see the manual page for class Integral.
 * @author Dmitry Chirkin
 * Ported to C++ by Jens Dreyer
 */

#ifndef FINDROOT_H_INCLUDED
#define FINDROOT_H_INCLUDED

#include <cmath>
#include "PROPOSAL/MathModel.h"
#include "PROPOSAL/DFunctionOfx.h"

class FindRoot:public MathModel
{
private:

    double _precision;
    int _maxSteps;
//----------------------------------------------------------------------------//

public:
    DFunctionOfx *function2use;

    //Constructors

	FindRoot(void);
//----------------------------------------------------------------------------//

	FindRoot(int maxSteps, double precision);

//----------------------------------------------------------------------------//
    //Destructor

	~FindRoot(void);
//----------------------------------------------------------------------------//
    //Memberfunctions

	double findRoot(double min, double max, double startX, DFunctionOfx *function2use, double rightSide);

//----------------------------------------------------------------------------//

	double function(double x);

//----------------------------------------------------------------------------//

	double dFunction(double x);
  


};

#endif // _FINDROOT_H_INCLUDED
