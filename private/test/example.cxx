
#include <fstream>
#include <iostream>

#include "PROPOSAL/PROPOSAL.h"

using namespace PROPOSAL;

int main(int argc, const char* argv[])
{
    int statistics = 100;

    std::cout << "c++ standard: " << __cplusplus << std::endl;

    if (argc >= 2)
    {
        int result;
        std::istringstream ss(argv[1]);
        ss.imbue(std::locale::classic());
        ss >> result;

        statistics = result;
        std::cout << "propagate " << statistics << " particles" << std::endl;
    }

    /**************************************************************************
     *              One Propagator for each particle definition                *
     **************************************************************************/

    Propagator prop_mu(MuMinusDef::Get(), "resources/config_ice.json");
    Propagator prop_tau(TauMinusDef::Get(), "resources/config_ice.json");

    // Therefor its needed to get the internal created particle first
    Particle& particle_mu  = prop_mu.GetParticle();
    Particle& particle_tau = prop_tau.GetParticle();

    particle_mu.SetEnergy(1e7); // [MeV]
    particle_mu.SetPropagatedDistance(0);
    particle_mu.SetPosition(Vector3D(0, 0, 100));
    particle_mu.SetDirection(Vector3D(0, 0, -1));

    particle_tau.SetEnergy(1e8); // [MeV]
    particle_tau.SetPropagatedDistance(0);
    particle_tau.SetPosition(Vector3D(0, 0, 0));
    particle_tau.SetDirection(Vector3D(0, 0, -1));

    std::vector<DynamicData*> sec_mu_direct = prop_mu.Propagate();

    /**************************************************************************
     *                        Using propagator service                         *
     **************************************************************************/

    // Possibility to register propagator in a service
    PropagatorService prop_service;
    prop_service.RegisterPropagator(prop_mu);
    prop_service.RegisterPropagator(prop_tau);

    // Define Particles to propagate
    Particle mu(MuMinusDef::Get());
    Particle tau(TauMinusDef::Get());

    // Set particle properties
    mu.SetEnergy(1e8); // [MeV]
    mu.SetPropagatedDistance(0);
    mu.SetPosition(Vector3D(0, 0, 0));
    mu.SetDirection(Vector3D(0, 0, -1));

    tau.SetEnergy(1e8); // [MeV]
    tau.SetPropagatedDistance(0);
    tau.SetPosition(Vector3D(0, 0, 0));
    tau.SetDirection(Vector3D(0, 0, -1));

    // Use service to propagate different particle
    std::vector<DynamicData*> sec_mu_service  = prop_service.Propagate(mu);
    std::vector<DynamicData*> sec_tau_service = prop_service.Propagate(tau);

    Output::getInstance().ClearSecondaryVector();
}
