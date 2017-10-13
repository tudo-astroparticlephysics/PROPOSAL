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
#pragma once

#ifndef ROOTFINDER_H
#define ROOTFINDER_H

#include <boost/function.hpp>

namespace PROPOSAL{

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
    bool operator==(const RootFinder &finder) const;
    bool operator!=(const RootFinder &finder) const;
    RootFinder(int maxSteps, double precision);

//----------------------------------------------------------------------------//
    //Memberfunctions

    /**
    * returns the value of the root bracketed between min and max.
    * Starting value of x is determined by 0&lt;=startX&lt;=1
    */

    double FindRoot(double min,
                    double max,
                    double startX,
                    boost::function<double (double)> function,
                    boost::function<double (double)> differentiated_function) const;

//----------------------------------------------------------------------------//

    void swap(RootFinder &finder);

//----------------------------------------------------------------------------//
    //Getter
	int GetMaxSteps() const {
		return maxSteps_;
	}
//----------------------------------------------------------------------------//
	double GetPrecision() const {
		return precision_;
	}
//----------------------------------------------------------------------------//
    //Setter
	void SetMaxSteps(int maxSteps);
	void SetPrecision(double precision);

//----------------------------------------------------------------------------//
    //Destructor
    ~RootFinder();
};

}

#endif // _ROOTFINDER_H
