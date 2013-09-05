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
    cout << *p << endl;
}
