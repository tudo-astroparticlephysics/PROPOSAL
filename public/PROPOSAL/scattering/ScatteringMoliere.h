#pragma once

#ifndef SCATTERING_MOLIERE_H
#define SCATTERING_MOLIERE_H

// #include <vector>
// #include <cmath>

#include "PROPOSAL/math/MathModel.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/PROPOSALParticle.h"
#include "PROPOSAL/scattering/Scattering.h"

namespace PROPOSAL {

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

    // Memberfunctions

    //----------------------------------------------------------------------------//

    // void SetBetaSq(double b) { betaSq_ = b; }

    // void SetChi0(std::vector<double> b) { chi0_ = b; }
    // void SetChi0(unsigned int i, double b) { chi0_.at(i) = b; }
    // void SetChiASq(std::vector<double> b) { chiASq_ = b; }
    // void SetChiASq(unsigned int i, double b) { chiASq_.at(i) = b; }
    // void SetChiCSq(double b) { chiCSq_ = b; }
    // void SetB(std::vector<double> b) { B_ = b; }
    // void SetB(unsigned int i, double b) { B_.at(i) = b; }

    //----------------------------------------------------------------------------//

    // double GetBetaSq() { return betaSq_; }

    // std::vector<double> GetWeigth() { return weight_; }

    // std::vector<double> GetChi0() { return chi0_; }
    // std::vector<double> GetChiASq() { return chiASq_; }
    // double GetChiCSq() { return chiCSq_; }
    // std::vector<double> GetB() { return B_; }

    // std::vector<double> GetChiCSqrtBq()
    // {
    //     std::vector<double> chiCBSq(numComp_);
    //     for (int i        = 0; i < numComp_; i++)
    //         chiCBSq.at(i) = sqrt(chiCSq_ * B_.at(i));

    //     return chiCBSq;
    // }

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
};
}

#endif // SCATTERING_MOLIERE_H
