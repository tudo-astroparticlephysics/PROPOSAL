#include <iostream>
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/Interpolant.h"
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

int main(){

    ifstream in;
    in.open("bin/Brems_dEdx.txt");

    char firstLine[256];
    in.getline(firstLine,256);
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


    while(in.good())
    {
        in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;
        cout<<para<<"\t"<<ecut<<"\t"<<vcut<<"\t"<<lpm<<"\t"<<energy<<"\t"<<med<<"\t"<<particleName<<"\t"<<dEdx<<"\t";

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);


        brems->SetParametrization(para);
        brems->EnableLpmEffect(lpm);
        brems->EnableDEdxInterpolation();
        if(para==1 && ecut==500 && vcut == -1 && lpm ==true){
            para = para*para;
        }
        dEdx_new=brems->CalculatedEdx();
        cout<<dEdx_new<<"\t"<<particle->GetMass()<<"\t"<<brems->GetVMax()<<"\t"<<brems->GetVMin()<<"\t"<<brems->GetVUp()<<endl;
        //ASSERT_NEAR(dEdx_new, dEdx, 5e-1*dEdx);

        delete cuts;
        delete medium;
        delete particle;
        delete brems;



    }

    return 0;
}

