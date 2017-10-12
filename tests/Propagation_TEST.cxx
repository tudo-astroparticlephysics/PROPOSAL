
// #include <iostream>
// #include <string>
// #include <cmath>
#include <vector>
#include <cmath>

#include "gtest/gtest.h"

#include "PROPOSAL/PROPOSAL.h"

using namespace std;
using namespace PROPOSAL;

TEST(Propagation , Test_nan) {

    int statistic = 10;
    int EmaxLog10 = 8;

    // Define Particles
    Particle mu(MuMinusDef::Get());
    Particle tau(TauMinusDef::Get());

    // One Propagator for each particle definition
    // medium/detector configuration
    Propagator prop_mu(mu.GetParticleDef(), "../src/resources/config_ice.json");
    Propagator prop_tau(tau.GetParticleDef(), "../src/resources/config_ice.json");

    // Possibility to register propagator in a service
    PropagatorService prop_service;
    prop_service.RegisterPropagator(prop_mu);
    prop_service.RegisterPropagator(prop_tau);

    std::vector<unsigned int> length_sec;

    for (int i = 0; i < statistic; i++)
    {
        // ----------------------------------------------------------------- //
        // Using propagator service
        // ----------------------------------------------------------------- //

        // Set particle properties
        mu.SetEnergy(pow(10,EmaxLog10));
        mu.SetPropagatedDistance(0);
        mu.SetPosition(Vector3D(0, 0, 0));
        mu.SetDirection(Vector3D(0, 0, -1));

        tau.SetEnergy(pow(10,EmaxLog10));
        tau.SetPropagatedDistance(0);
        tau.SetPosition(Vector3D(0, 0, 0));
        tau.SetDirection(Vector3D(0, 0, -1));

        // Use service to propagate different particle
        std::vector<DynamicData*> sec_mu = prop_service.Propagate(mu);
        std::vector<DynamicData*> sec_tau = prop_service.Propagate(tau);

        // ----------------------------------------------------------------- //
        // Using propagator directly
        // ----------------------------------------------------------------- //

        // Therefor its needed to get the internal created particle first
        Particle& particle = prop_mu.GetParticle();

        particle.SetEnergy(pow(10,EmaxLog10));
        particle.SetPropagatedDistance(0);
        particle.SetPosition(Vector3D(0, 0, 0));
        particle.SetDirection(Vector3D(0, 0, -1));

        std::vector<DynamicData*> sec_mu_direct = prop_mu.Propagate();

        length_sec.push_back(sec_mu_direct.size());
    }

    for(std::vector<unsigned int>::iterator it = length_sec.begin(); it != length_sec.end(); ++it)
    {
        std::cout << "length of secondies: " << *it << std::endl;
    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}




