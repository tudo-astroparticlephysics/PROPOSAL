#include <iostream>
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/Interpolant.h"
#include "PROPOSAL/Ionization.h"

using namespace std;

double integrate(double min, double max, int N, double (*func)(double) ){
    double dx   =   max-min;
    dx          /=  N;
    double  integral    = func(min) + func(max);
            integral    /=2;

    for(int i=1;i<N;i++)integral += func(dx*i+min);

    return integral*dx;
}

double X2(double x){
    return x*x;
}

double X_YY(double x, double y){
    return x + y*y;
}

double _3X_2(double x){
    return 3*x + 2;
}


double dNdx_dNdxRnd(double rnd, double energy){

    double dNdxRnd;
    double dNdx;

    Medium *medium = new Medium("antares water",1.);
    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
    particle->SetEnergy(energy);
    EnergyCutSettings *cuts = new EnergyCutSettings(500.,0.1);

    CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);

    brems->SetParametrization(1.);

    dNdx=brems->CalculatedNdx();
    dNdxRnd=brems->CalculatedNdx(rnd);
    cout<<"dNdx\t"<<dNdx<<"\t"<<"dNdxRnd\t "<<dNdxRnd<<"\t";

    return dNdx-dNdxRnd;

}


int main(){

    double energy1;
    double rnd;
    double diff;


    char firstLine[256];

    double dEdx_new;
    double energy;
    double dEdx;
    double ecut;
    double vcut;
    string med;
    string particleName;
    bool lpm;
    int para;

    cout.precision(16);

    Medium *medium = new Medium("ice",1.);
    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
    particle->SetEnergy(10000);
    EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);

    CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);
    brems->SetParametrization(1);
    brems->EnableLpmEffect(0);

    double rnd1= 0.9655968558508899;
    double rnd2= 0.0596232083626091;

    double gna = brems->CalculateStochasticLoss(rnd1,rnd2);
    brems->EnableDNdxInterpolation();
    double gna2 = brems->CalculateStochasticLoss(rnd1,rnd2);

    cout << "should be: " << 9262.664586574878 << endl;
    cout << "brems(" << rnd1 << "," << rnd2 << "): " << gna << endl;
    cout << "Interpol: " << gna2 << endl;
    cout<<gna/gna2<<endl;



    return 0;
}

