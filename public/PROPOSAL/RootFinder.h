/**
 * Original comment of the Java class by Dmitry Chirkin:
 * This class provides routines for calculating roots by
 * the combination of the Newton-Raphson method and bisection.
 * Methods contained here are based on the Numerical Recipes
 * (W. H. Press et al.).
 *
 * @author Dmitry Chirkin
 * Ported to C++ Jan-Hendrik Köhne
 */

#ifndef ROOTFINDER_H
#define ROOTFINDER_H

#include "boost/function.hpp"

class RootFinder
{
private:

    int     maxSteps_;
    double  precision_;


public:

    //Constructors

    RootFinder();
    RootFinder(const RootFinder&);
    RootFinder& operator=(const RootFinder&);
    RootFinder(int maxSteps, double precision);

//----------------------------------------------------------------------------//
    //Destructor

    ~RootFinder();
//----------------------------------------------------------------------------//
    //Memberfunctions

    double FindRoot(double min,
                    double max,
                    double startX,
                    boost::function<double (double)> function,
                    boost::function<double (double)> differentiated_function,
                    double rightSide);

	int GetMaxSteps() const {
		return maxSteps_;
	}

	double GetPrecision() const {
		return precision_;
	}

	void SetMaxSteps(int maxSteps);
	void SetPrecision(double precision);
};

#endif // _ROOTFINDER_H
