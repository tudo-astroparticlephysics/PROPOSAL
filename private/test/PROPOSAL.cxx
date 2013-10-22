#include <iostream>
#include <fstream>
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/Interpolant.h"
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/Epairproduction.h"
#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/ContinuousRandomization.h"
#include "PROPOSAL/Geometry.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/StandardNormal.h"

//Stuff for LOG4CPLUS
#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/configurator.h>
#include <log4cplus/layout.h>
#include <iomanip>

using namespace log4cplus;
using namespace std;



int main(int argc, char** argv)
{
/*
    Propagator *pr = new Propagator("resources/configuration");
    pr->set_seed(1234);
    cout << "blah!" << endl;

    Particle *p = new Particle("mu");
    p->SetEnergy(9e3);
    pr->SetParticle(p);

    pr->Propagate(p);


    const int init1[] = {1000, 1, 0, 0};
     const int init2[] = {1000, 1, 0, 0};
     const int init3[] = {1000, 1, 0, 1};
     std::vector<int> vector1 ( init1, init1 + 4 );
     std::vector<int> vector2 ( init2, init2 + 4 );
     std::vector<int> vector3 ( init3, init3 + 4 );
     bool is_equal = false;
     if ( vector1.size() < vector2.size() )
       is_equal = std::equal ( vector1.begin(), vector1.end(), vector2.begin() );
     else
       is_equal = std::equal ( vector2.begin(), vector2.end(), vector1.begin() );
     std::cout<<"operator==: "<< ( vector1 == vector2 ) <<'\n';
     std::cout<<"std::equal: "<< is_equal <<'\n';
          std::cout<<"operator==: "<< ( vector1 == vector3 ) <<'\n';
     vector3 = vector2;
     std::cout<<"operator==: "<< ( vector1 == vector3 ) <<'\n';
*/


    // Test ScatteringMoliere

    // scattering angle

    ScatteringMoliere *scat1 = new ScatteringMoliere();

    Particle *part1 = new Particle("mu",0,0,0,0,0,7.3e3,0);
    Medium *medi1 = new Medium("salt", 1.);

    scat1->Scatter(1.44,part1,medi1);

    ofstream out1;
    out1.open("test_salt_getrndm.txt");

    double divisor = 1;
    for(int i = 0; i < 1e7/divisor; i++)
    {
        out1 << scat1->GetRandom() << std::endl;
    }

    out1.close();

    // now : theta

    ScatteringMoliere *scat2 = new ScatteringMoliere();

    Medium *medi2 = new Medium("salt", 1.);
    Particle *part2;


    ofstream out2;
    out2.open("test_salt_theta.txt");

    for(int i = 0; i < 1e6/divisor; i++)
    {
        part2 = new Particle("mu",0,0,0,0,0,7.3e3,0);

        scat2->Scatter(1.44,part2,medi2);

        out2 << part2->GetTheta() << std::endl;
    }

    out2.close();



    // Test ScatteringFirstOrder (Highland)
    ScatteringFirstOrder *scat = new ScatteringFirstOrder();

    Particle *part = new Particle("mu",0,0,0,0,0,7.3e3,0);
    Medium *medi = new Medium("salt", 1.);

    std::cout << scat->CalculateTheta0(1.44,part, medi ) << std::endl;


    ofstream out;
    out.open("test_firstorder.txt");

    double X;

    double x = 1./0.;
    double Zero=0.;
    double Nan = sqrt(-1.);
    cerr << "x: " << x << endl;
    cerr << "x/0. : " << x/Zero << endl;
    cerr << "x*0 : " << x*Zero << endl;
    cerr << "Sqrt(-1) : " << Nan << endl;

    for(int i = 0; i < 1e6/divisor; i++)
    {
        part = new Particle("mu",0,0,0,0,0,7.3e3,0);

        scat->Scatter(1.44,part, medi );

        X = part->GetTheta();

        if(X == X) out << X << std::endl;
    }

    out.close();
}


