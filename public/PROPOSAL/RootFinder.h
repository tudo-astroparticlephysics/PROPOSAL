
/******************************************************************************
 *																			  *
 * This file is part of the simulation tool PROPOSAL.						  *
 *																			  *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,		  *
 * 				      Chair Experimental Physics 5b							  *
 *																			  *
 * This software may be modified and distributed under the terms of a		  *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE". 									  *
 *																			  *
 * Modifcations to the LGPL License:										  *
 *																			  *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the		  *
 *         following reference:												  *
 *																			  *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001										  *
 *																			  *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *		   GitHub webpage													  *
 *																			  *
 *		   "https://github.com/tudo-astroparticlephysics/PROPOSAL"			  *
 *																			  *
 ******************************************************************************/

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
                    boost::function<double (double)> differentiated_function,
                    double rightSide);

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
