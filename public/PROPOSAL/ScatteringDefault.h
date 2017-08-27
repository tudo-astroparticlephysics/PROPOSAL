
#pragma once

#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Interpolant.h"
#include "PROPOSAL/Scattering.h"

namespace PROPOSAL {

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

    void Scatter(PROPOSALParticle&, const std::vector<CrossSections*>&, double dr, double ei, double ef);
    long double CalculateTheta0(const PROPOSALParticle&,
                                const std::vector<CrossSections*>&,
                                double dr,
                                double ei,
                                double ef);
    void EnableInterpolation(const PROPOSALParticle&, const std::vector<CrossSections*>&, std::string path = "");
    void DisableInterpolation();

    private:
    ScatteringDefault& operator=(const ScatteringDefault&); // Undefined & not allowed

    double FunctionToIntegral(const PROPOSALParticle&, const std::vector<CrossSections*>&, double energy);
    double FunctionToBuildInterpolant(const PROPOSALParticle&, const std::vector<CrossSections*>&, double energy);

    bool do_interpolation_;
    int order_of_interpolation_;

    Integral integral_;
    Interpolant* interpolant_;
    Interpolant* interpolant_diff_;
};

}
