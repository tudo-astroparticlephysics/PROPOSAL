
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

#ifndef SCATTERING_H
#define SCATTERING_H

// #include <vector>
// #include <string>

#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/Interpolant.h"
#include "PROPOSAL/Integral.h"
// #include "PROPOSAL/MathModel.h"
// #include "PROPOSAL/PROPOSALParticle.h"

namespace PROPOSAL{

/**
  * \brief This class provides the scattering routine provided by moliere.
  *
  * More precise scattering angles will be added soon.
  */


class Scattering : public MathModel
{


private:

    double x0_;

    bool do_interpolation_;
    int order_of_interpolation_;

    Integral* integral_;
    Interpolant* interpolant_;
    Interpolant* interpolant_diff_;

    PROPOSALParticle* particle_;
    PROPOSALParticle* backup_particle_;
    std::vector<CrossSections*> crosssections_;

//----------------------------------------------------------------------------//

public:

    /**
     * \brief Default Constructor
     *
     * Constructor which sets "default" settings.
     */
    Scattering();

//----------------------------------------------------------------------------//

    Scattering(std::vector<CrossSections*> crosssections);
//----------------------------------------------------------------------------//

    Scattering(const Scattering&);
    Scattering& operator=(const Scattering&);
    bool operator==(const Scattering &scattering) const;
    bool operator!=(const Scattering &scattering) const;
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
    // Memberfunctions

    long    double  CalculateTheta0(double dr, double ei, double ef);
    void            Scatter(double dr, double ei, double ef);
    double          FunctionToIntegral(double energy);
    double          FunctionToBuildInterpolant(double energy);
    void            EnableInterpolation(std::string path = "");
    void            DisableInterpolation();
//----------------------------------------------------------------------------//

    void swap(Scattering &scattering);

//----------------------------------------------------------------------------//

    //Setter

    void SetParticle(PROPOSALParticle* particle);
    void SetCrosssections(std::vector<CrossSections*> crosssections);
//----------------------------------------------------------------------------//
    // Getter

    PROPOSALParticle* GetParticle()
    {
        return particle_;
    }

    double GetX0()
    {
        return x0_;
    }

//----------------------------------------------------------------------------//
    // destructors
    ~Scattering() {}

    PROPOSALParticle *GetBackup_particle() const;
    void SetBackup_particle(PROPOSALParticle *backup_particle);
    void RestoreBackup_particle();

};

}


#endif //SCATTERING_H
