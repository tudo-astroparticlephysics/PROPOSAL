
#pragma once

#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/scattering/Scattering.h"

namespace PROPOSAL {

class CrossSection;
/**
  * \brief This class provides the scattering routine provided by moliere.
  *
  * More precise scattering angles will be added soon.
  */
class ScatteringDefault: public Scattering
{
    public:
    ScatteringDefault();
    ScatteringDefault(const ScatteringDefault&);
    ~ScatteringDefault() {}

    virtual Scattering* clone() const { return new ScatteringDefault(*this); }
    static Scattering* create() { return new ScatteringDefault(); }

    // bool operator==(const ScatteringDefault& scattering) const;
    // bool operator!=(const ScatteringDefault& scattering) const;
    // void swap(ScatteringDefault& scattering);

    private:
    ScatteringDefault& operator=(const ScatteringDefault&); // Undefined & not allowed

    RandomAngles CalculateRandomAngle(const PROPOSALParticle&, const Medium&, double dr, double disp);
    long double CalculateTheta0(const PROPOSALParticle&,
                                const Medium&,
                                double dr,
                                double disp);

};

}
