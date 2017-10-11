#pragma once

#ifndef SCATTERING_MOLIERE_H
#define SCATTERING_MOLIERE_H

// #include <vector>
// #include <cmath>

// #include "PROPOSAL/medium/Medium.h"
// #include "PROPOSAL/particle/PROPOSALParticle.h"
#include "PROPOSAL/scattering/Scattering.h"

namespace PROPOSAL {

// class PROPOSALParticle;
// class Medium;

class ScatteringMoliere : public Scattering
{
    public:
    // constructor
    ScatteringMoliere(PROPOSALParticle&, const Medium&);
    ScatteringMoliere(PROPOSALParticle&, const ScatteringMoliere&);
    ScatteringMoliere(const ScatteringMoliere&);
    ~ScatteringMoliere();

    Scattering* clone() const { return new ScatteringMoliere(*this); }
    virtual Scattering* clone(PROPOSALParticle& particle) const { return new ScatteringMoliere(particle, *this); }
    static Scattering* create(PROPOSALParticle& particle, const Medium& medium) { return new ScatteringMoliere(particle, medium); }

    // ScatteringMoliere& operator=(const ScatteringMoliere&);
    // bool operator==(const ScatteringMoliere& scattering) const;
    // bool operator!=(const ScatteringMoliere& scattering) const;
    //----------------------------------------------------------------------------//

    //----------------------------------------------------------------------------//
    //----------------------------------------------------------------------------//

    // destructors

    private:
    ScatteringMoliere& operator=(const ScatteringMoliere&); // Undefined & not allowed

    RandomAngles CalculateRandomAngle(double dr, double ei, double ef);

    const Medium* medium_;

    // double dx_;     // traversing thickness in cm
    // double betaSq_; // beta² = v²/c²
    // double p_;      // momentum in MeV/c
    // double m_;      // mass in MeV/c²

    // medium
    // Medium* medium_;
    int numComp_;                // number of components in medium
    std::vector<double> Zi_;     // nuclear charge of different components
    std::vector<double> ki_;     // number of atoms in molecule of different components
    std::vector<double> Ai_;     // atomic number of different components
    double A_;
    std::vector<double> weight_; // mass weights of different components

    // scattering parameters
    // std::vector<double> chi0_;
    // std::vector<double> chiASq_; // screening angle² in rad²
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
