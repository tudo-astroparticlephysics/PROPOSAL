/*! \file   Scattering.h
*   \brief  Header file for the Scattering routines.
*
*   For more details see the class documentation.
*
*   \date   2013.06.13
*   \author Tomasz Fuchs
*/
#pragma once

// #include <vector>
// #include <string>

#include "PROPOSAL/Scattering.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/Medium.h"

namespace PROPOSAL{

/**
  * \brief This class provides the scattering routine provided by moliere.
  *
  * More precise scattering angles will be added soon.
  */
class ScatteringFirstOrder : public Scattering
{
    public:
    ScatteringFirstOrder();
    ScatteringFirstOrder(const ScatteringFirstOrder&);
    ~ScatteringFirstOrder();

    virtual Scattering* clone() const { return new ScatteringFirstOrder(*this); }

    void Scatter(PROPOSALParticle&, const std::vector<CrossSections*>&, double dr, double ei, double ef);
    // Do nothing, not interpolation for scattering moliere
    virtual void EnableInterpolation(const PROPOSALParticle&, const std::vector<CrossSections*>&, std::string path = "");
    virtual void DisableInterpolation();

    private:
    ScatteringFirstOrder& operator=(const ScatteringFirstOrder&); // Undefined & not allowed

    double CalculateTheta0(const PROPOSALParticle&, const Medium&, double dr);
};

}
