#ifndef SCATTERING_MOLIERE_H
#define SCATTERING_MOLIERE_H

#include "vector"
#include <cmath>

#include "PROPOSAL/Medium.h"
#include "PROPOSAL/Particle.h"

#include "PROPOSAL/MathModel.h"
#include "PROPOSAL/Integral.h"

using namespace std;

class ScatteringMoliere
{
private:

    //particle
    double dx;                  //traversing thickness in cm
    double betaSq;              //beta² = v²/c²
    double p;                   //momentum in MeV/c
    double m;                   //mass in MeV/c²

    //medium
    Medium* medium;
    int numComp;                //number of components in medium
    vector<double> Z;           //nuclear charge of different components
    vector<double> ki;          //number of atoms in molecule of different components
    vector<double> Ai;          //atomic number of different components
    vector<double> weight;      //mass weights of different components


    //scattering parameters
    vector<double> chi0;
    vector<double> chiASq;      //screening angle² in rad²
    double chiCSq;              //characteristic angle² in rad²
    vector<double> B;

    double thetaMax;              //maximum randomly generated angle


    MathModel* MathMachine;
    Integral* IntMachine;


public:

    //constructor
    ScatteringMoliere();

    //Old constructor
    //ScatteringMoliere(double dr, Particle* part, Medium* med);

//----------------------------------------------------------------------------//

    // Memberfunctions
    void Scatter(double dr, Particle* part, Medium* med);

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

    void CalcBetaSq();

    void CalcWeight();

    void CalcChi0();
    void CalcChiASq();
    void CalcChiCSq();
    void CalcB();

//----------------------------------------------------------------------------//

    void SetBetaSq(double b) { betaSq = b; }

    void SetChi0(vector<double> b) { chi0 = b; }
    void SetChi0(unsigned int i, double b) { chi0.at(i) = b; }
    void SetChiASq(vector<double> b) { chiASq = b; }
    void SetChiASq(unsigned int i, double b) { chiASq.at(i) = b; }
    void SetChiCSq(double b) { chiCSq = b; }
    void SetB(vector<double> b) { B = b; }
    void SetB(unsigned int i, double b) { B.at(i) = b; }

    void SetThetaMax(double b) { thetaMax = b; }

//----------------------------------------------------------------------------//

    double GetBetaSq() { return betaSq; }

    vector<double> GetWeigth() { return weight; }

    vector<double> GetChi0() { return chi0; }
    vector<double> GetChiASq() { return chiASq; }
    double GetChiCSq() { return chiCSq; }
    vector<double> GetB() { return B; }

    vector<double> GetChiCSqrtBq()
    {
        vector<double> chiCBSq(medium->GetNumComponents());
        for(int i = 0; i < medium->GetNumComponents() ; i++) chiCBSq.at(i) = sqrt(chiCSq*B.at(i));

        return chiCBSq;
    }

    double GetThetaMax() { return thetaMax; }

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

    double f1M(double x);
    double f2M(double x);

    double f(double theta);

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

    int BinarySearch(int n, const double *array, double value);
    double GetRandom();

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

    // destructors
    ~ScatteringMoliere() {}
};

#endif // SCATTERING_MOLIERE_H
