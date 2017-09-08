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
class CrossSection;

class Scattering
{
    public:
    Scattering();
    virtual ~Scattering();

    virtual Scattering* clone() const = 0; // virtual constructor idiom (used for deep copies)

    void Scatter(PROPOSALParticle&, const std::vector<CrossSection*>&, double dr, double ei, double ef);

    virtual void EnableInterpolation(const PROPOSALParticle&, const std::vector<CrossSection*>&, std::string path = "") = 0;
    virtual void DisableInterpolation() = 0;

    protected:
    struct RandomAngles
    {
        double sx, sy, tx, ty;
    };

    virtual RandomAngles CalculateRandomAngle(const PROPOSALParticle&, const std::vector<CrossSection*>&, double dr, double ei, double ef) = 0;
};

}
