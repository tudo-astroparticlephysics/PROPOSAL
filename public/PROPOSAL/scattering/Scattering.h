/*! \file   Scattering.h
*   \brief  Header file for the Scattering bug routines.
*
*   This version has a major bug and produces too small scattering angles.
*
*   \date   2013.08.19
*   \author Tomasz Fuchs
*/
#pragma once

#include <vector>
#include <string>

namespace PROPOSAL {

class PROPOSALParticle;
class Medium;

class Scattering
{
    public:
    Scattering();
    virtual ~Scattering();

    virtual Scattering* clone() const = 0; // virtual constructor idiom (used for deep copies)

    void Scatter(PROPOSALParticle&, const Medium&, double dr, double disp);

    protected:
    struct RandomAngles
    {
        double sx, sy, tx, ty;
    };

    virtual RandomAngles CalculateRandomAngle(const PROPOSALParticle&, const Medium&, double dr, double disp) = 0;
};

}
