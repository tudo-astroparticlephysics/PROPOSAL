/*! \file   Scattering.h
*   \brief  Header file for the Scattering routines.
*
*   For more details see the class documentation.
*
*   \date   
*   \author 
*/
#pragma once

#include "PROPOSAL/scattering/Scattering.h"

namespace PROPOSAL{

/**
  * \brief This class provides the scattering routine provided by moliere.
  *
  * More precise scattering angles will be added soon.
  */
class ScatteringNoScattering : public Scattering
{
    public:
    ScatteringNoScattering(PROPOSALParticle&, const Medium&);
    ScatteringNoScattering(PROPOSALParticle&, const ScatteringNoScattering&);
    ScatteringNoScattering(const ScatteringNoScattering&);
    ~ScatteringNoScattering();

    virtual Scattering* clone() const { return new ScatteringNoScattering(*this); }
    virtual Scattering* clone(PROPOSALParticle& particle) const { return new ScatteringNoScattering(particle, *this); }
    static Scattering* create(PROPOSALParticle& particle, const Medium& medium) { return new ScatteringNoScattering(particle, medium); }

    private:
    ScatteringNoScattering& operator=(const ScatteringNoScattering&); // Undefined & not allowed

    RandomAngles CalculateRandomAngle(double dr, double ei, double ef);

    const Medium* medium_;
};

}
