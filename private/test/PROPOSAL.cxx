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
    cout.precision(16);
    cout << "baslsadkadsadasdh" << endl;
    double min, max;
    min = 1;
    max = 3;
    int N = 10000000;
    cout << "min: " << min << endl;
    cout << "max: " << max << endl;
    cout << "--------------------------" << endl;
    cout << "integral x*x: "        << integrate(min,max,N,X2)      << endl;
    cout << "integral 3*x + 2: "    << integrate(min,max,N,_3X_2)   << endl;

    Medium *med = new Medium("ice",1.);
    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
    Bremsstrahlung *brems = new Bremsstrahlung();
    brems->SetMedium(med);
    brems->EnableLpmEffect(true);

    for(int i = 3 ; i < 13 ; i++){
        brems->SetParametrization(1);
        particle->SetEnergy(pow(10.,i));
        brems->SetParticle(particle);
        cout<<brems->ElasticBremsstrahlungCrossSection(0.2,0)<<endl;
    }

    ifstream in;
    in.open("Brems_dEdx.txt");
    string aps;
    in>>aps;
    cout<<aps<<endl;
    brems->SetParametrizationLimit(1.);

    Integral* Int = new Integral();
//    cout << "IntegralKlasse x*x: " << Int->IntegrateClosed(0,3,X2) << endl;
    bool testInterpolant = true;
    if(testInterpolant){
        int max = 100;
        double xmin = 0;
        double xmax = 20;
        int romberg = 5;
        bool rational = false;
        bool relative = false;
        bool isLog = false;
        int rombergY = 5;
        bool rationalY = false;
        bool relativeY = false;
        bool logSubst = false;
        Interpolant* Pol1 = new Interpolant(max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);
        double SearchX = 3;
        cout << "Interpolating f(" << SearchX << "): " << Pol1->interpolate(3) << endl;


                int max2 = 100;
        double x2min = 0;
        double x2max = 20;
        int romberg2 = 5;
        bool rational2 = false;
        bool relative2 = false;
        bool isLog2 = false;

                Interpolant* Pol2 = new Interpolant(max,  xmin, xmax, max2, x2min, x2max, X_YY,
                                                    romberg, rational, relative, isLog,
                                                    romberg2, rational2, relative2, isLog2,
                                                    rombergY, rationalY, relativeY, logSubst);

                cout << "x+y*y (2,3): " <<  Pol2->interpolate(2.,3.) << endl;
    }


    return 0;
}

