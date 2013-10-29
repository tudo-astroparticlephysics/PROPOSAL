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
    Medium *med = new Medium("ice",1);
        EnergyCutSettings *ecuts = new EnergyCutSettings();

        Particle *part = new Particle("mu");
        part->SetEnergy(10e6);
        part->SetPropagatedDistance(5e5);

        Propagator *prop = new Propagator("resources/configuration_IceOnly");


        prop->set_seed(142);
        prop->SetParticle(part);


        prop->Propagate(part);

        std::cout << "GetX = " << part->GetX() << std::endl;
        std::cout << "GetY = " << part->GetY() << std::endl;
        std::cout << "GetZ = " << part->GetZ() << std::endl;

}


