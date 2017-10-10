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

#include "PROPOSAL/scattering/Scattering.h"
#include "PROPOSAL/particle/PROPOSALParticle.h"
#include "PROPOSAL/medium/Medium.h"

namespace PROPOSAL{

/**
  * \brief This class provides the scattering routine provided by moliere.
  *
  * More precise scattering angles will be added soon.
  */
class ScatteringFirstOrder : public Scattering
{
    public:
    ScatteringFirstOrder(PROPOSALParticle&, const Medium&);
    ScatteringFirstOrder(PROPOSALParticle&, const ScatteringFirstOrder&);
    ScatteringFirstOrder(const ScatteringFirstOrder&);
    ~ScatteringFirstOrder();

    virtual Scattering* clone() const { return new ScatteringFirstOrder(*this); }
    virtual Scattering* clone(PROPOSALParticle& particle) const { return new ScatteringFirstOrder(particle, *this); }
    static Scattering* create(PROPOSALParticle& particle, const Medium& medium) { return new ScatteringFirstOrder(particle, medium); }

    private:
    ScatteringFirstOrder& operator=(const ScatteringFirstOrder&); // Undefined & not allowed

    RandomAngles CalculateRandomAngle(double dr, double ei, double ef);
    double CalculateTheta0(double dr);

    const Medium* medium_;
};

}
