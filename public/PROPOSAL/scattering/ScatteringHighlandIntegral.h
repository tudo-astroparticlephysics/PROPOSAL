
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
class ScatteringHighlandIntegral : public Scattering
{
public:
    ScatteringHighlandIntegral(Particle&, const Utility&);
    ScatteringHighlandIntegral(Particle&, const Utility&, const InterpolationDef&);

    // Copy constructor
    ScatteringHighlandIntegral(Particle&, const Utility&, const ScatteringHighlandIntegral&);
    ScatteringHighlandIntegral(const ScatteringHighlandIntegral&);
    ~ScatteringHighlandIntegral();

    virtual Scattering* clone() const { return new ScatteringHighlandIntegral(*this); }
    virtual Scattering* clone(Particle& particle, const Utility& utility) const
    {
        return new ScatteringHighlandIntegral(particle, utility, *this);
    }

private:
    ScatteringHighlandIntegral& operator=(const ScatteringHighlandIntegral&); // Undefined & not allowed

    bool compare(const Scattering&) const;

    RandomAngles CalculateRandomAngle(double dr, double ei, double ef);
    long double CalculateTheta0(double dr, double ei, double ef);

    UtilityDecorator* scatter_;
};

} // namespace PROPOSAL
