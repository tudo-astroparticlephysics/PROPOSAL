
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
    ScatteringDefault(Particle&, Utility&);
    ScatteringDefault(Particle&, Utility&, InterpolationDef);
    ScatteringDefault(Particle&, const ScatteringDefault&);

    ScatteringDefault(const ScatteringDefault&);
    ~ScatteringDefault();

    virtual Scattering* clone() const { return new ScatteringDefault(*this); }
    virtual Scattering* clone(Particle& particle) const { return new ScatteringDefault(particle, *this); }
    static Scattering* create(Particle& particle, Utility& utility) { return new ScatteringDefault(particle, utility); }

    // bool operator==(const ScatteringDefault& scattering) const;
    // bool operator!=(const ScatteringDefault& scattering) const;
    // void swap(ScatteringDefault& scattering);

    private:
    ScatteringDefault& operator=(const ScatteringDefault&); // Undefined & not allowed

    RandomAngles CalculateRandomAngle(double dr, double ei, double ef);
    long double CalculateTheta0(double dr,
                                double ei,
                                double ef);

    UtilityDecorator* scatter;
};

}
