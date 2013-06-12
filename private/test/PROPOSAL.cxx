#include <iostream>
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/Interpolant.h"
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/Epairproduction.h"
#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/ContinuousRandomization.h"
#include "PROPOSAL/Geometry.h"

using namespace std;

//double integrate(double min, double max, int N, double (*func)(double) ){
//    double dx   =   max-min;
//    dx          /=  N;
//    double  integral    = func(min) + func(max);
//            integral    /=2;

//    for(int i=1;i<N;i++)integral += func(dx*i+min);

//    return integral*dx;
//}

//double X2(double x){
//    return x*x;
//}

//double X_YY(double x, double y){
//    return x + y*y;
//}

//double _3X_2(double x){
//    return 3*x + 2;
//}


//double dNdx_dNdxRnd(double rnd, double energy){

//    double dNdxRnd;
//    double dNdx;

//    Medium *medium = new Medium("antares water",1.);
//    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
//    particle->SetEnergy(energy);
//    EnergyCutSettings *cuts = new EnergyCutSettings(500.,0.1);

//    CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);

//    brems->SetParametrization(1.);

//    dNdx=brems->CalculatedNdx();
//    dNdxRnd=brems->CalculatedNdx(rnd);
//    cout<<"dNdx\t"<<dNdx<<"\t"<<"dNdxRnd\t "<<dNdxRnd<<"\t";

//    return dNdx-dNdxRnd;

//}

//double loss(double energy,double rnd1,double rnd2)
//{
//    Medium *medium = new Medium("antares water",1.);
//    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
//    particle->SetEnergy(energy);
//    EnergyCutSettings *cuts = new EnergyCutSettings(500.,0.1);

//    CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);

//    brems->SetParametrization(1.);


//    return brems->CalculateStochasticLoss(rnd1,rnd2);
//}


//int main(){


//    char firstLine[256];

//    double dEdx_new;
//    double dEdx;
//    double ecut = 500;
//    double vcut = 0.05;
//    string med = "uranium";
//    string particleName = "mu";
//    bool lpm = 0;
//    int para;

//    cout.precision(16);


//    particleName = "mu";
//    ecut = 500;
//    vcut = 0.001;
//    med = "ice";
//    lpm = 0;
//    double energy = 1E4;

//    Medium *medium = new Medium(med,1.);
//    Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
//    particle->SetEnergy(energy);
//    EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);


////    CrossSections *epair= new Epairproduction(particle, medium, cuts);
////    epair->EnableLpmEffect(lpm);

//    double rnd1= 0.9655968558508899;
//    double rnd2= 0.0596232083626091;

////    cout << "should be: " << 2192.996821337525 << endl;
////    double gna = epair->CalculateStochasticLoss(rnd1,rnd2);
////    cout << "e: " << gna << endl;

////    cout << "Interpolating..." << endl;
////    epair->EnableDNdxInterpolation();
////    cout << "should be: " << 2190.419541590142 << endl;
////    double gna2 = epair->CalculateStochasticLoss(rnd1,rnd2);
////    cout << "Interpol: " << gna2 << endl;
////    cout<<gna/gna2<<endl;

//    double neu;
//    double alt;


//    //Bremsstrahlung *brems= new Bremsstrahlung(particle, medium, cuts);
//    Ionization *brems = new Ionization(particle, medium, cuts);


//    for(int k =0; k<1000;k++)
//    {
//    for( int i =1; i<11;i++)
//    {
//         for(int j =1; j<11;j++)
//         {
//            energy=120 * pow(10,k/100.);
//            brems->GetParticle()->SetEnergy(energy);
//            rnd1=1./i;
//            rnd2 =1./j;

//            //brems->CalculatedNdx(rnd1);
//            neu=brems->CalculateStochasticLossNew(rnd1, rnd2);

//            alt=brems->CalculateStochasticLoss(rnd1, rnd2);
//            if(alt==alt && fabs((alt-neu)/alt)>1e-5)
//            cout<<rnd1<<"\t"<<rnd2<<"--\t "<<alt<<"\t "<<neu<<"\t"<<(alt-neu)/alt<<"\t--"<<endl;
//            //cout<<alt<<endl;
//            //cout<<neu<<endl;

//         }
//    }
//    }




//    return 0;
//}

int main(int argc, char** argv){

    double x,y,z,theta,phi;
    x = -15.8385;
    y = -8.7836;
    z = -3.99414;
    theta = 157.057;
    phi = 317.619;
    Particle * particle = new Particle("mu",x,y,z,theta,phi,0,0);
    Geometry A;

    double x0,y0,z0,radius,inner_radius,height;
    x0=-4.23397;
    y0=-5.88271;
    z0=-14.7617;
    radius = 10;
    inner_radius = 6.72594;
    height = 10;
    A.InitCylinder(x0,y0,z0,radius,inner_radius,height);
    A.IsParticleInside(particle);

}
