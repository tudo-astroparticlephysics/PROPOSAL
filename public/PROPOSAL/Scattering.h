/*! \file   Scattering.h
*   \brief  Header file for the Scattering bug routines.
*
*   This version has a major bug and produces too small scattering angles.
*
*   \date   2013.08.19
*   \author Tomasz Fuchs
*/
#pragma once

#ifndef SCATTERING_H
#define SCATTERING_H

// #include <vector>
// #include <string>

#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/Interpolant.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/MathModel.h"
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
    bool do_interpolation_;
    int order_of_interpolation_;

    Integral* integral_;
    Interpolant* interpolant_;
    Interpolant* interpolant_diff_;

public:

    /**
     * \brief Default Constructor
     *
     * Constructor which sets "default" settings.
     */
    Scattering();

//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//

    Scattering(const Scattering&);
    Scattering& operator=(const Scattering&);
    bool operator==(const Scattering &scattering) const;
    bool operator!=(const Scattering &scattering) const;
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
    // Memberfunctions

    void            Scatter(PROPOSALParticle&, const std::vector<CrossSections*>&, double dr, double ei, double ef);

    long    double  CalculateTheta0(const PROPOSALParticle&, const std::vector<CrossSections*>&, double dr, double ei, double ef);
    double          FunctionToIntegral(const PROPOSALParticle&, const std::vector<CrossSections*>&, double energy);
    double          FunctionToBuildInterpolant(const PROPOSALParticle&, const std::vector<CrossSections*>&, double energy);

    void            EnableInterpolation(const PROPOSALParticle&, const std::vector<CrossSections*>&, std::string path = "");
    void            DisableInterpolation();

//----------------------------------------------------------------------------//

    void swap(Scattering &scattering);

//----------------------------------------------------------------------------//
    // destructors
    ~Scattering() {}
};

}


#endif //SCATTERING_H
