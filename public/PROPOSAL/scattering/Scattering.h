/*! \file   Scattering.h
*   \brief  Header file for the Scattering bug routines.
*
*   This version has a major bug and produces too small scattering angles.
*
*   \date   2013.08.19
*   \author Tomasz Fuchs
*/
#pragma once

namespace PROPOSAL {

class Particle;
class Utility;

class Scattering
{
    public:
    Scattering(Particle&);
    Scattering(const Scattering&);
    virtual ~Scattering();

    bool operator==(const Scattering& scattering) const;
    bool operator!=(const Scattering& scattering) const;

    virtual Scattering* clone() const = 0; // virtual constructor idiom (used for deep copies)
    virtual Scattering* clone(Particle&, const Utility&) const = 0; // virtual constructor idiom (used for deep copies)

    void Scatter(double dr, double ei, double ef);

    const Particle& GetParticle() const { return particle_; }

    protected:
    Scattering& operator=(const Scattering&); // Undefined & not allowed

    // Implemented in child classes to be able to use equality operator
    virtual bool compare(const Scattering&) const = 0;

    struct RandomAngles
    {
        double sx, sy, tx, ty;
    };

    virtual RandomAngles CalculateRandomAngle(double dr, double ei, double ef) = 0;

    Particle& particle_;
};

}
