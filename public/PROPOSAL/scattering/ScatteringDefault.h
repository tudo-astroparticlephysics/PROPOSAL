
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

    void EnableInterpolation(const PROPOSALParticle&, const std::vector<CrossSection*>&, std::string path = "");
    void DisableInterpolation();


    private:
    RandomAngles CalculateRandomAngle(const PROPOSALParticle&, const std::vector<CrossSection*>&, double dr, double ei, double ef);
    long double CalculateTheta0(const PROPOSALParticle&,
                                const std::vector<CrossSection*>&,
                                double dr,
                                double ei,
                                double ef);
    ScatteringDefault& operator=(const ScatteringDefault&); // Undefined & not allowed

    double FunctionToIntegral(const PROPOSALParticle&, const std::vector<CrossSection*>&, double energy);
    double FunctionToBuildInterpolant(const PROPOSALParticle&, const std::vector<CrossSection*>&, double energy);

    bool do_interpolation_;
    int order_of_interpolation_;

    Integral integral_;
    Interpolant* interpolant_;
    Interpolant* interpolant_diff_;
};

}
