
#include <iostream>

#include "PROPOSAL/PROPOSAL.h"

using namespace PROPOSAL;

int main()
{
    // Define Particles
    ParticleDef particle_def = MuMinusDef::Get();
    particle_def.mass = 1000; // [MeV]
    PROPOSALParticle mu(particle_def);
    PROPOSALParticle tau(TauMinusDef::Get());

    // One Propagator for each particle definition
    // medium/detector configuration
    Propagator prop_mu(mu.GetParticleDef(), "config.json");
    Propagator prop_tau(tau.GetParticleDef(), "config.json");

    // Possibility to register propagator in a service
    PropagatorService prop_service;
    prop_service.RegisterPropagator(prop_mu);
    prop_service.RegisterPropagator(prop_tau);

    // Using propagator service

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
    std::vector<DynamicData*> sec_mu = prop_service.Propagate(mu);
    std::vector<DynamicData*> sec_tau = prop_service.Propagate(tau);

    // Using propagator directly

    // Therefor its needed to get the internal created particle first
    PROPOSALParticle& particle = prop_mu.GetParticle();

    particle.SetEnergy(1e8); // [MeV]
    particle.SetPropagatedDistance(0);
    particle.SetPosition(Vector3D(0, 0, 0));
    particle.SetDirection(Vector3D(0, 0, -1));

    std::vector<DynamicData*> sec_mu_direct = prop_mu.Propagate();

    // for(std::vector<unsigned int>::iterator it = length_sec.begin(); it != length_sec.end(); ++it)
    // {
    //     std::cout << "length of secondies: " << *it << std::endl;
    // }
}
