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

//    double x,y,z,theta,phi;
//    x = -15.8385;
//    y = -8.7836;
//    z = -3.99414;
//    theta = 157.057;
//    phi = 317.619;
//    Particle * particle = new Particle("mu",x,y,z,theta,phi,0,0);
//    Geometry A;

//    double x0,y0,z0,radius,inner_radius,height;
//    x0=-4.23397;
//    y0=-5.88271;
//    z0=-14.7617;
//    radius = 10;
//    inner_radius = 6.72594;
//    height = 10;
//    A.InitCylinder(x0,y0,z0,radius,inner_radius,height);
//    A.IsParticleInside(particle);

//    Particle *mu = new Particle("mu");
//    Medium* med = new Medium("air",1);
//    EnergyCutSettings* cut = new EnergyCutSettings(100,-1);
//    Bremsstrahlung* br1 = new Bremsstrahlung(mu,med,cut);
//    br1->EnableDNdxInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaaah");
//    br1->EnableDEdxInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaaah");

//    Bremsstrahlung* br2 = new Bremsstrahlung(mu,med,cut);
//    br2->EnableDNdxInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaaah");
//    br2->EnableDEdxInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaah");

//    Particle *mu = new Particle("mu");
//    Medium* med = new Medium("air",1);
//    EnergyCutSettings* cut = new EnergyCutSettings(100,-1);
//    Ionization* io1 = new Ionization(mu,med,cut);
//    io1->EnableDNdxInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaaah");
//    io1->EnableDEdxInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaaah");

//    Ionization* io2 = new Ionization(mu,med,cut);
//    io2->EnableDNdxInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaaah");
//    io2->EnableDEdxInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaaah");

//    Particle *mu = new Particle("mu");
//    Medium* med = new Medium("air",1);
//    EnergyCutSettings* cut = new EnergyCutSettings(100,-1);
//    Epairproduction* ep1 = new Epairproduction(mu,med,cut);
//    ep1->EnableDNdxInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaaah");
//    ep1->EnableDEdxInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaaah");

//    Epairproduction* ep2 = new Epairproduction(mu,med,cut);
//    ep2->EnableDNdxInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaaah");
//    ep2->EnableDEdxInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaaah");



//    Particle *mu = new Particle("mu");
//    Medium* med = new Medium("air",1);
//    EnergyCutSettings* cut = new EnergyCutSettings(100,0.001);
//    Photonuclear* p1 = new Photonuclear(mu,med,cut);
//    p1->EnableDNdxInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaaah");
//    p1->EnableDEdxInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaaah");

//    Photonuclear* p2 = new Photonuclear(mu,med,cut);
//    p2->EnableDNdxInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaaah");
//    p2->EnableDEdxInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaaah");

//    Particle *mu = new Particle("mu");
//    Medium* med = new Medium("air",1);
//    EnergyCutSettings* cut = new EnergyCutSettings(100,0.001);
//    ProcessCollection* p1 = new ProcessCollection(mu,med,cut);
//    p1->EnableContinuousRandomization();
//    p1->EnableInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaaah");

//    ProcessCollection* p2 = new ProcessCollection(mu,med,cut);
//    p2->EnableContinuousRandomization();
//    p2->EnableInterpolation("/home/koehne/PROPOSAL_restructure/buildrestructure_Mar_11_2013/uaaaaah");

//    Geometry *k = new Geometry();
//    k->InitBox(0,0,2,3,5,6);
//    cout <<*k<<endl;

//    EnergyCutSettings *cut = new EnergyCutSettings();
//    cout<<*cut<<endl;

//    Medium *med = new Medium("antares_water",1.003);
//    cout<<*med<<endl;

//    Particle *p = new Particle(12,13,"mu",1,2,3,30,43,1e6,0.6,0);
//    cout<<*p<<endl;
    Propagator *pr = new Propagator("resources/configuration");
//    Propagator *pr1 = new Propagator();
//    pr1->EnableInterpolation("resources/tables");
    double x;

    ofstream out;
    out.open("dist_mu_9e6.txt");
    int number_of_particles =   1e6;
    for(int i =0 ;i< number_of_particles; i++)
    {
        if(i%10000==0)cout<<1.*i/1e4<<endl;
        pr->set_seed(i);
        Particle *p = new Particle("mu");
        p->SetEnergy(9e6);
        pr->SetParticle(p);
        x   =   pr->Propagate(p);
out<<x<<endl;
    }


//    Propagator* prop = new Propagator();

//    prop->EnableInterpolation("/data/LocalApps/LocalFiles/tables");
//    prop->GetParticle()->SetEnergy(702);
//    double distance = prop->Propagate(1E10);


//    double startenergy = 370;
//    Medium *medium = new Medium("ice",1.);
//    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
//    particle->SetEnergy(startenergy);
//    EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);

//    std::vector<CrossSections*> vecOfProcColl;
//    CrossSections* ion = new Ionization(particle,medium,cuts);

//    CrossSections* brems = new Bremsstrahlung(particle,medium,cuts);
//    brems->SetParametrization(1);

//    CrossSections* epair = new Epairproduction(particle,medium,cuts);

//    CrossSections* photo = new Photonuclear(particle,medium,cuts);
//    photo->SetParametrization(12);


//    vecOfProcColl.push_back(ion);
//    vecOfProcColl.push_back(brems);
//    vecOfProcColl.push_back(epair);
//    vecOfProcColl.push_back(photo);

//    Scattering* scat = new Scattering(vecOfProcColl);
//    scat->GetParticle()->SetPhi(0);
//    scat->GetParticle()->SetTheta(0);
//    scat->GetParticle()->SetX(0);
//    scat->GetParticle()->SetY(0);
//    scat->GetParticle()->SetZ(0);

//    scat->EnableInterpolation("/data/LocalApps/LocalFiles/tables");

//    scat->Scatter(50,startenergy,MMU);

//    cout << *(scat->GetParticle()) << endl;
//    }

}
