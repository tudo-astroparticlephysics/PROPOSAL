#pragma once

#ifndef SCATTERING_MOLIERE_H
#define SCATTERING_MOLIERE_H

#include <vector>

#include "PROPOSAL/scattering/Scattering.h"

namespace PROPOSAL {

class Medium;

class ScatteringMoliere : public Scattering
{
public:
    // constructor
    ScatteringMoliere(Particle&, const Medium&);
    ScatteringMoliere(Particle&, const ScatteringMoliere&);
    ScatteringMoliere(const ScatteringMoliere&);
    ~ScatteringMoliere();

    Scattering* clone() const { return new ScatteringMoliere(*this); }
    virtual Scattering* clone(Particle& particle, const Utility& utility) const
    {
        (void)utility;
        return new ScatteringMoliere(particle, *this);
    }

private:
    ScatteringMoliere& operator=(const ScatteringMoliere&); // Undefined & not allowed

    bool compare(const Scattering&) const;

    RandomAngles CalculateRandomAngle(double dr, double ei, double ef);

    const Medium* medium_;

    int numComp_; // number of components in medium
    double ZSq_A_average_;
    std::vector<double> Zi_;        // nuclear charge of different components
    std::vector<double> weight_ZZ_; // mass weights of different components time Z^2
    double weight_ZZ_sum_;          // inverse of sum of mass weights of different components time Z^2
    int max_weight_index_;          // index of the maximium of mass weights of different components

    // scattering parameters
    double chiCSq_; // characteristic angle² in rad²
    std::vector<double> B_;

    //----------------------------------------------------------------------------//
    //----------------------------------------------------------------------------//

    double f1M(double x);
    double f2M(double x);

    double f(double theta);

    double F1M(double x);
    double F2M(double x);

    double F(double theta);

    //----------------------------------------------------------------------------//
    //----------------------------------------------------------------------------//

    double GetRandom(double pre_factor);
};
} // namespace PROPOSAL

#endif // SCATTERING_MOLIERE_H
