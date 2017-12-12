
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

#ifndef MATHMODEL_H_
#define MATHMODEL_H_

// #include <fstream>

#include <boost/function.hpp>
#include <boost/random.hpp>


namespace PROPOSAL{

class MathModel
{

public:

//----------------------------------------------------------------------------//
    //Constructor

    MathModel();
    MathModel(const MathModel&);
    MathModel& operator=(const MathModel&);

//----------------------------------------------------------------------------//
    //Destructor

    virtual ~MathModel();

//----------------------------------------------------------------------------//
    //Memberfunctions

    double RandomDouble();

//----------------------------------------------------------------------------//

    static void set_seed(int seed);

    /** @brief Set a custom random number generator
     *
     * Classes that contain other subclasses of MathModel should
     * override this to pass the new RNG on to their members.
     */
    virtual void SetRandomNumberGenerator(boost::function<double ()> &f);

private:
        static double DefaultRandomDouble();

    boost::function<double ()> rng_;
    static boost::mt19937 *default_rng_;
};

}

#endif // MATHMODEL_H_
