

/*
 * Interpolate.h
 *
 *  Created on: 02.08.2010
 *      Author: koehne
 */

#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

#include <vector>
#include <cmath>

#include "MathModel.h"
#include "FunctionInt2.h"
#include "FunctionInt.h"


#include "Output.h"


/**
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

class Interpolate: public MathModel,  public FunctionInt{ /// implements FunctionInt{

private:
	int romberg, rombergY;
        double *iX;
        double *iY;
        double *c;
        double *d;
	int max;
	double xmin, xmax, step;
	bool rational, relative;

	FunctionInt2 *function2int;
        Interpolate* interpolate_;
        int row, starti;
	bool rationalY, relativeY;

	bool reverse, self, flag;		// Self is setted to true in constructor
	bool isLog, logSubst;

	double precision, worstX;
	double precision2, worstX2;
	double precisionY, worstY;

	bool fast;			// Is setted to true in constructor			

	double x_save, y_save;	// Is setted to 1 and 0 in constructor

	Output *output;

        /**
         * interpolates f(x) based on the values iY[i]=f(iX[i]) in the romberg-vicinity of x
         */

            double interpolate(double x, int start);

        //----------------------------------------------------------------------------------------------------//



protected:

	const static double bigNumber=-300;
        const static double aBigNumber=-299;


    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes class for the 1-dimensional function
     */

public:

	//----------------------------------------------------------------------------------------------------//

	/* Default Constructor */

	Interpolate();


    Interpolate(int max, double xmin, double xmax, FunctionInt *function2int,
		       int romberg, bool rational, bool relative, bool isLog,
		       int rombergY, bool rationalY, bool relativeY, bool logSubst);

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes class for the 2-dimensional function
     */

	Interpolate(int max1, double x1min, double x1max, int max2, double x2min, double x2max, FunctionInt2 *function2int,
		       int romberg1, bool rational1, bool relative1, bool isLog1,
		       int romberg2, bool rational2, bool relative2, bool isLog2,
		       int rombergY, bool rationalY, bool relativeY, bool logSubst);

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes class for the 1-dimensional function if the arrays already exist.
     */

        Interpolate(std::vector<double> x, std::vector<double> y, int romberg, bool rational, bool relative);

    //----------------------------------------------------------------------------------------------------//

    /**
     * auxiliary class initializer
     */

       void InitInterpolate(int max, double xmin, double xmax,
			int romberg, bool rational, bool relative, bool isLog,
			int rombergY, bool rationalY, bool relativeY, bool logSubst);

    //----------------------------------------------------------------------------------------------------//

    /**
     * defines a function for every row
     */

	double functionInt(double x);

    //----------------------------------------------------------------------------------------------------//

    /**<<"iX\t"
     * interpolates f(x) for 1d function
     */

	double interpolate(double x);

    //----------------------------------------------------------------------------------------------------//

    /**
     * interpolates f(x) for 2d function
     */

	double interpolate(double x1, double x2);

    //----------------------------------------------------------------------------------------------------//

    /**
     * interpolates f(x) for 1d function if the arrays already exist.
     */

	double interpolateArray(double x);

    //----------------------------------------------------------------------------------------------------//



    /**
     * finds x: f(x)=y, 1d initialization required
     */

	double findLimit(double y);

    //----------------------------------------------------------------------------------------------------//

    /**
     * finds x: f(a,x)=y, 2d initialization required
     */

	double findLimit(double x1, double y);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Exp if not zero
     */

	double exp(double x);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Log if not zero
     */

	double log(double x);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Log it again
     */

	double slog(double x);

    //----------------------------------------------------------------------------------------------------//

        Interpolate* get_interpolate(){return interpolate_;}

        //std::vector<double> get_iX(){return iX;}
        //std::vector<double> get_iY(){return iY;}

        int get_max(){return max;}
};

#endif /* INTERPOLATE_H_ */
