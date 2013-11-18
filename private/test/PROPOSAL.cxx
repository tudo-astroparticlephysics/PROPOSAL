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
//    Output::getInstance().EnableROOTOutput("test_Output.root");
    Output::getInstance().EnableASCIIOutput("ascii_test");
    for(int i = 0; i< (int)(1e4) ; i++)
    {
        Particle* part = new Particle(i,i,"mu",0,0,0,0,0,0,0,0);
        part->SetEnergy(1e6);

        propa->Propagate(part);
    }

    Output::getInstance().Close();





/*
    double alpha,blobel_alpha;

    alpha = 0.68;
    blobel_alpha = alpha/2 + 0.5;

    blobel_alpha = 0.84;
    alpha = (blobel_alpha -0.5)*2;

    cerr << "alpha:\t" << alpha  << endl;
    cerr << "B_alpha:\t" << blobel_alpha << endl;

    double Z = SQRT2*erfInv(alpha);
    cerr << "Z(" << alpha << "): " << Z << endl;

    double low, high;
    double event = 1;
    double epsilonH = 0.46;
    double epsilonL = 0.48;

    for(event = 1; event < 11; event++)
    {
        high = pow( Z/2 + sqrt(event + 0.5 + epsilonH),2.);
        low  = pow(Z /2 - sqrt(event + 0.5 - epsilonL), 2.);

        cerr << "k: \t" << event << "\t\tLimits: [" << low << " - " << high << "]\n";
    }
*/
}


