// FindRoot header file


/**
 * This class provides routines for calculating roots by the combination of the Newton-Raphson method and bisection.
 * Include the function to be integrated in a class that implements the interface FunctionOfx (defined below).
 * Methods contained here are based on the Numerical Recipes (W. H. Press et al.).
 * <pre>
 * interface DFunctionOfx extends FunctionOfx{
 *     double dFunction(double x);
 * }
 * </pre>
 * For the definition of interface FunctionOfx see the manual page for class Integral.
 * @author Dmitry Chirkin
 */

#ifndef FINDROOT_H_INCLUDED
#define FINDROOT_H_INCLUDED

#include <cmath>
#include "MathModel.h"
#include "DFunctionOfx.h"

class FindRoot:public MathModel{
public:
	FindRoot(void);
	FindRoot(int maxSteps, double precision);
	~FindRoot(void);
	DFunctionOfx *function2use;

	double findRoot(double min, double max, double startX, DFunctionOfx *function2use, double rightSide);
	double function(double x);
	double dFunction(double x);
  

private:
	double _precision;
    	int _maxSteps;
    };

#endif // _FINDROOT_H_INCLUDED
