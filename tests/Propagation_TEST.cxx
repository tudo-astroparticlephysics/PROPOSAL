
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

    // one propagator for each particle definition
    // medium/detector configuration
    Propagator prop_mu(mu.GetParticleDef(), "resources/config_ice.json");
    Propagator prop_tau(tau.GetParticleDef(), "resources/config_ice.json");

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

TEST(Propagation , particle_type)
{
    ifstream in;
    string filename = "bin/TestFiles/propagation.txt";
    in.open(filename.c_str());

    if (!in.good())
    {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    // Just skip the header
    char firstLine[256];
    in.getline(firstLine,256);

    int statistic = 100;
    double energy = 1e8;

    in >> statistic >> energy;
    std::cout << statistic << std::endl;
    std::cout << energy << std::endl;

    Propagator prop_mu(MuMinusDef::Get(), "resources/config_ice.json");
    Particle& mu = prop_mu.GetParticle();

    mu.SetEnergy(energy);
    mu.SetPropagatedDistance(0);
    mu.SetPosition(Vector3D(0, 0, 0));
    mu.SetDirection(Vector3D(0, 0, -1));

    std::vector<std::string> names;
    std::vector<double> lengths;
    std::vector<double> energies;
    std::vector<Vector3D> positions;
    std::vector<Vector3D> directions;

    RandomGenerator::Get().SetSeed(0);

    // Buffer to read data in

    std::string name;
    double length = 0.0;
    double sec_energy = 0.0, x = 0.0, y = 0.0, z = 0.0, dx = 0.0, dy = 0.0, dz = 0.0;

    double error = 1e-5;

    // Read

    cout.precision(16);

    for (int i = 0; i < statistic; ++i)
    {
        mu.SetEnergy(energy);
        mu.SetPropagatedDistance(0);
        mu.SetPosition(Vector3D(0, 0, 0));
        mu.SetDirection(Vector3D(0, 0, -1));

        std::vector<DynamicData*> sec_mu_direct = prop_mu.Propagate();

        std::string name_new = "";

        for (unsigned int j = 0; j < sec_mu_direct.size(); ++j)
        {
            if (sec_mu_direct[j]->GetTypeId() == DynamicData::Particle)
            {
                Particle* particle = dynamic_cast<Particle*>(sec_mu_direct[j]);
                // names.push_back(particle->GetName());
                name_new = particle->GetName();
            }
            else
            {
                switch (sec_mu_direct[j]->GetTypeId())
                {
                    case DynamicData::Brems:
                        // names.push_back("Brems");
                        name_new = "Brems";
                        break;
                    case DynamicData::Epair:
                        // names.push_back("Epair");
                        name_new = "Epair";
                        break;
                    case DynamicData::NuclInt:
                        // names.push_back("NuclInt");
                        name_new = "NuclInt";
                        break;
                    case DynamicData::DeltaE:
                        // names.push_back("DeltaE");
                        name_new = "DeltaE";
                        break;
                    default:
                        break;
                }
            }

            in >> name >> length >> sec_energy >> x >> y >> z >> dx >> dy >> dz;

            double energy_new = sec_mu_direct[j]->GetEnergy();
            double lenght_new = sec_mu_direct[j]->GetPropagatedDistance();
            Vector3D position = sec_mu_direct[j]->GetPosition();
            Vector3D direction = sec_mu_direct[j]->GetDirection();

            EXPECT_TRUE(name == name_new);
            ASSERT_NEAR(length, lenght_new, error*length);
            ASSERT_NEAR(sec_energy, energy_new, error * sec_energy);
            ASSERT_NEAR(x, position.GetX(), error * std::abs(x));
            ASSERT_NEAR(y, position.GetY(), error * std::abs(y));
            ASSERT_NEAR(z, position.GetZ(), error * std::abs(z));
            ASSERT_NEAR(dx, direction.GetX(), error * std::abs(dx));
            ASSERT_NEAR(dy, direction.GetY(), error * std::abs(dy));
            ASSERT_NEAR(dz, direction.GetZ(), error * std::abs(dz));
        }
    }

    in.close();
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}




