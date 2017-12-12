
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

#ifndef SCATTERING_FIRSTORDER_H
#define SCATTERING_FIRSTORDER_H

// #include <vector>
// #include <string>

#include "PROPOSAL/MathModel.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/Medium.h"

namespace PROPOSAL{

/**
  * \brief This class provides the scattering routine provided by moliere.
  *
  * More precise scattering angles will be added soon.
  */


class ScatteringFirstOrder : public MathModel
{

//----------------------------------------------------------------------------//

public:

    /**
     * \brief Default Constructor
     *
     * Constructor which sets "default" settings.
     */
    ScatteringFirstOrder();


//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
    // Memberfunctions

    double  CalculateTheta0(double dr, PROPOSALParticle* part, Medium* med);
    void    Scatter(double dr, PROPOSALParticle* part, Medium* med);


//----------------------------------------------------------------------------//
    // destructors
    ~ScatteringFirstOrder() {}


};

}


#endif //SCATTERING_FIRSTORDER_H
