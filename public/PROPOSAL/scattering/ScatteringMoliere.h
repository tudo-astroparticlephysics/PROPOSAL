#pragma once

#ifndef SCATTERING_MOLIERE_H
#define SCATTERING_MOLIERE_H

// #include <vector>
// #include <cmath>

// #include "PROPOSAL/medium/Medium.h"
// #include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/scattering/Scattering.h"

namespace PROPOSAL {

// class Particle;
// class Medium;

class ScatteringMoliere : public Scattering
{
    public:
    // constructor
    ScatteringMoliere(Particle&, const Medium&);
    ScatteringMoliere(Particle&, const ScatteringMoliere&);
    ScatteringMoliere(const ScatteringMoliere&);
    ~ScatteringMoliere();

    Scattering* clone() const { return new ScatteringMoliere(*this); }
    virtual Scattering* clone(Particle& particle) const { return new ScatteringMoliere(particle, *this); }
    static Scattering* create(Particle& particle, const Medium& medium) { return new ScatteringMoliere(particle, medium); }

    // ScatteringMoliere& operator=(const ScatteringMoliere&);
    // bool operator==(const ScatteringMoliere& scattering) const;
    // bool operator!=(const ScatteringMoliere& scattering) const;
    //----------------------------------------------------------------------------//

    //----------------------------------------------------------------------------//
    //----------------------------------------------------------------------------//

    // destructors

    private:
    ScatteringMoliere& operator=(const ScatteringMoliere&); // Undefined & not allowed

    bool compare(const Scattering&) const;

    RandomAngles CalculateRandomAngle(double dr, double ei, double ef);
    // dr is the traversing thickness in cm

    const Medium* medium_;

    int numComp_;                // number of components in medium
    double A_average_;
    std::vector<double> Zi_;     // nuclear charge of different components
    std::vector<double> weight_; // mass weights of different components

    // scattering parameters
    double chiCSq_;              // characteristic angle² in rad²
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

    double GetRandom();

    //----------------------------------------------------------------------------//
    //-----------------------------Coefficients-----------------------------------//
    //---for calculating the power series approximation of the moliere function---//

    // static const double c1[100];
    // static const double c2[100];
    // static const double c2large[50];
    // static const double s2large[50];
    // static const double C1large[50];

};
}

#endif // SCATTERING_MOLIERE_H
