
#include <iostream>

#include "PROPOSAL/PROPOSAL.h"

using namespace PROPOSAL;

int main()
{
    // Define Particles
    // ParticleDef particle_def = MuMinusDef::Get();
    // particle_def.mass = 1000; // [MeV]
    // Particle mu(MuMinusDef::Get());
    Particle tau(TauMinusDef::Get());

    // One Propagator for each particle definition
    // medium/detector configuration

    // Propagator prop_mu(mu.GetParticleDef(), "resources/config_ice.json");
    Propagator prop_tau(tau.GetParticleDef(), "resources/config_ice.json");

    // Possibility to register propagator in a service
    // PropagatorService prop_service;
    // prop_service.RegisterPropagator(prop_mu);
    // prop_service.RegisterPropagator(prop_tau);

    // Using propagator service

    // Set particle properties
    // mu.SetEnergy(1e8); // [MeV]
    // mu.SetPropagatedDistance(0);
    // mu.SetPosition(Vector3D(0, 0, 0));
    // mu.SetDirection(Vector3D(0, 0, -1));

    // tau.SetEnergy(1e8); // [MeV]
    // tau.SetPropagatedDistance(0);
    // tau.SetPosition(Vector3D(0, 0, 0));
    // tau.SetDirection(Vector3D(0, 0, -1));

    // Use service to propagate different particle
    // std::vector<DynamicData*> sec_mu = prop_service.Propagate(mu);
    // std::vector<DynamicData*> sec_tau = prop_service.Propagate(tau);

    // Using propagator directly

    // Therefor its needed to get the internal created particle first
    // Particle& particle_mu = prop_mu.GetParticle();
    Particle& particle_tau = prop_tau.GetParticle();

    // particle_mu.SetEnergy(1e8); // [MeV]
    // particle_mu.SetPropagatedDistance(0);
    // particle_mu.SetPosition(Vector3D(0, 0, 0));
    // particle_mu.SetDirection(Vector3D(0, 0, -1));

    particle_tau.SetEnergy(1e8); // [MeV]
    particle_tau.SetPropagatedDistance(0);
    particle_tau.SetPosition(Vector3D(0, 0, 0));
    particle_tau.SetDirection(Vector3D(0, 0, -1));

    // std::cout << "propagate muon" << std::endl;
    // std::vector<DynamicData*> sec_mu_direct = prop_mu.Propagate();
    // std::cout << sec_mu_direct.size() << std::endl;

    std::cout << "propagate tau" << std::endl;
    std::vector<DynamicData*> sec_tau_direct = prop_tau.Propagate();
    std::cout << sec_tau_direct.size() << std::endl;

}
