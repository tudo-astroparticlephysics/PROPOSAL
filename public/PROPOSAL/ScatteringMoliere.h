#ifndef SCATTERING_MOLIERE_H
#define SCATTERING_MOLIERE_H

#include "vector"
#include <cmath>

#include "PROPOSAL/Medium.h"
#include "PROPOSAL/Particle.h"

#include "PROPOSAL/MathModel.h"

using namespace std;

class ScatteringMoliere
{
private:

    //particle
    double dx_;                  //traversing thickness in cm
    double betaSq_;              //beta² = v²/c²
    double p_;                   //momentum in MeV/c
    double m_;                   //mass in MeV/c²

    //medium
    Medium* medium_;
    int numComp_;                //number of components in medium
    vector<double> Zi_;          //nuclear charge of different components
    vector<double> ki_;          //number of atoms in molecule of different components
    vector<double> Ai_;          //atomic number of different components
    vector<double> weight_;      //mass weights of different components


    //scattering parameters
    vector<double> chi0_;
    vector<double> chiASq_;      //screening angle² in rad²
    double chiCSq_;              //characteristic angle² in rad²
    vector<double> B_;


    MathModel* MathMachine_;


public:

    //constructor
    ScatteringMoliere();
    ScatteringMoliere(const ScatteringMoliere &);
    ScatteringMoliere& operator=(const ScatteringMoliere&);
    bool operator==(const ScatteringMoliere &scattering) const;
    bool operator!=(const ScatteringMoliere &scattering) const;
//----------------------------------------------------------------------------//

    // Memberfunctions
    void Scatter(double dr, Particle* part, Medium* med);
    void swap(ScatteringMoliere &scattering);
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

    void CalcBetaSq();

    void CalcWeight();

    void CalcChi0();
    void CalcChiASq();
    void CalcChiCSq();
    void CalcB();

//----------------------------------------------------------------------------//

    void SetBetaSq(double b) { betaSq_ = b; }

    void SetChi0(vector<double> b) { chi0_ = b; }
    void SetChi0(unsigned int i, double b) { chi0_.at(i) = b; }
    void SetChiASq(vector<double> b) { chiASq_ = b; }
    void SetChiASq(unsigned int i, double b) { chiASq_.at(i) = b; }
    void SetChiCSq(double b) { chiCSq_ = b; }
    void SetB(vector<double> b) { B_ = b; }
    void SetB(unsigned int i, double b) { B_.at(i) = b; }

//----------------------------------------------------------------------------//

    double GetBetaSq() { return betaSq_; }

    vector<double> GetWeigth() { return weight_; }

    vector<double> GetChi0() { return chi0_; }
    vector<double> GetChiASq() { return chiASq_; }
    double GetChiCSq() { return chiCSq_; }
    vector<double> GetB() { return B_; }

    vector<double> GetChiCSqrtBq()
    {
        vector<double> chiCBSq(medium_->GetNumComponents());
        for(int i = 0; i < medium_->GetNumComponents() ; i++) chiCBSq.at(i) = sqrt(chiCSq_*B_.at(i));

        return chiCBSq;
    }

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
//----------------------------------------------------------------------------//

    // destructors
    ~ScatteringMoliere() {}
};

#endif // SCATTERING_MOLIERE_H
