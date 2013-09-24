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

}
