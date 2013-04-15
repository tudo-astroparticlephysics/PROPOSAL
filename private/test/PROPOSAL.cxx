#include <iostream>
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/Interpolant.h"
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/Epairproduction.h"

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


    char firstLine[256];

    double dEdx_new;
    double dEdx;
    double ecut = 500;
    double vcut = -1;
    string med = "uranium";
    string particleName = "mu";
    bool lpm = 0;
    int para;

    cout.precision(16);


    particleName = "mu";
    ecut = 500;
    vcut = -1;
    med = "ice";
    lpm = 0;
    double energy = 1E4;

    Medium *medium = new Medium(med,1.);
    Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
    particle->SetEnergy(energy);
    EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

    CrossSections *epair= new Epairproduction(particle, medium, cuts);
    epair->EnableLpmEffect(lpm);

    double rnd1= 0.9655968558508899;
    double rnd2= 0.0596232083626091;

    cout << "should be: " << 5.105934174938547e-07 << endl;

    double gna = epair->CalculatedNdx(rnd1);
    cout << "dNdx(): " << gna << endl;

    epair->EnableDNdxInterpolation();
    double gna2 = epair->CalculatedNdx(rnd1);
    cout << "Interpol: " << gna2 << endl;
    cout<<gna/gna2<<endl;







    return 0;
}

