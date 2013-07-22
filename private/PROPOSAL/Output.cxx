#include "PROPOSAL/Output.h"

using namespace std;

vector<Particle*> Output::secondarys_;

void Output::SetLoggingConfigurationFile(std::string file)
{
    PropertyConfigurator::doConfigure(LOG4CPLUS_TEXT(file));
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Output::FillSecondaryVector(Particle *particle, int secondary_id, pair<double,string> energy_loss, double distance)
{
    string secondary_name;

    if(energy_loss.second.compare("Epairproduction")==0)
    {
        secondary_name  =    "epair";
    }
    else if(energy_loss.second.compare("Ionization")==0)
    {
        secondary_name  =    "delta";
    }
    else if(energy_loss.second.compare("Bremsstrahlung")==0)
    {
        secondary_name  =    "brems";
    }
    else if(energy_loss.second.compare("Photonuclear")==0)
    {
        secondary_name  =    "munu";
    }
    else  //decay
    {
        secondary_name  =   energy_loss.second;
    }

    Particle *particle_to_store   =   new Particle(particle->GetParentParticleId(), secondary_id, secondary_name, particle->GetX(), particle->GetY(), particle->GetZ(), particle->GetTheta(), particle->GetPhi(), energy_loss.first, particle->GetT(), distance);
    secondarys_.push_back(particle_to_store);

}
