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

class Medium;

/**
  * \brief This class provides the scattering routine provided by moliere.
  *
  * More precise scattering angles will be added soon.
  */
class ScatteringNoScattering : public Scattering
{
    public:
    ScatteringNoScattering(Particle&, const Medium&);
    ScatteringNoScattering(Particle&, const ScatteringNoScattering&);
    ScatteringNoScattering(const ScatteringNoScattering&);
    ~ScatteringNoScattering();

    virtual Scattering* clone() const { return new ScatteringNoScattering(*this); }
    virtual Scattering* clone(Particle& particle, const Utility& utility) const { (void) utility; return new ScatteringNoScattering(particle, *this); }

    private:
    ScatteringNoScattering& operator=(const ScatteringNoScattering&); // Undefined & not allowed

    bool compare(const Scattering&) const;

    RandomAngles CalculateRandomAngle(double dr, double ei, double ef);

    const Medium* medium_;
};

}
