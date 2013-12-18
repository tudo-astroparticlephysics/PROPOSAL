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

#include <time.h>
#include <boost/math/special_functions/erf.hpp>
#define erfInv(x)   boost::math::erf_inv(x)

using namespace log4cplus;
using namespace std;




int main(int argc, char** argv)
{

    Propagator* propa = new Propagator("resources/configuration");

    Output::getInstance().EnableROOTOutput("TestOutput.root");
    for(int i = 0; i< (int)(1e4) ; i++)
    {
        Particle* part = new Particle(i,i,"mu",0,0,0,0,0,0,0,0);
        part->SetEnergy(1e6);

        propa->Propagate(part);
    }

    Output::getInstance().Close();
}


