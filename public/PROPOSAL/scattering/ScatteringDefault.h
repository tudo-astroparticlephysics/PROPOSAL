
#pragma once

#include "PROPOSAL/scattering/Scattering.h"

namespace PROPOSAL {

class Utility;
class UtilityDecorator;
struct InterpolationDef;
/**
  * \brief This class provides the scattering routine provided by moliere.
  *
  * More precise scattering angles will be added soon.
  */
class ScatteringDefault: public Scattering
{
    public:
    ScatteringDefault(Particle&, const Utility&);
    ScatteringDefault(Particle&, const Utility&, const InterpolationDef&);

    // Copy constructor
    ScatteringDefault(Particle&, const Utility&, const ScatteringDefault&);
    ScatteringDefault(const ScatteringDefault&);
    ~ScatteringDefault();

    virtual Scattering* clone() const { return new ScatteringDefault(*this); }
    virtual Scattering* clone(Particle& particle, const Utility& utility) const { return new ScatteringDefault(particle, utility, *this); }

    // void swap(ScatteringDefault& scattering);

    private:
    ScatteringDefault& operator=(const ScatteringDefault&); // Undefined & not allowed

    bool compare(const Scattering&) const;

    RandomAngles CalculateRandomAngle(double dr, double ei, double ef);
    long double CalculateTheta0(double dr, double ei, double ef);

    UtilityDecorator* scatter_;
};

}
