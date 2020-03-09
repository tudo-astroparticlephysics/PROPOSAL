
#include <cmath>
#include <vector>

#include "gtest/gtest.h"

#include "PROPOSAL/PROPOSAL.h"

using namespace PROPOSAL;

TEST(Comparison, Comparison_equal)
{
    Sector::Definition sector_def;
    sector_def.location = Sector::ParticleLocation::InsideDetector;
    sector_def.SetMedium( std::make_shared<const Water>() );
    sector_def.SetGeometry(Sphere().create());
    sector_def.scattering_model            = ScatteringFactory::Moliere;
    sector_def.cut_settings                = EnergyCutSettings();
    sector_def.do_continuous_randomization = true;

    std::vector<Sector::Definition> sec_defs;
    sec_defs.push_back(sector_def);

    Propagator prop_a(MuMinusDef::Get(), sec_defs, Sphere().create());
    Propagator prop_b(MuMinusDef::Get(), sec_defs, Sphere().create());

    EXPECT_TRUE(prop_a == prop_b);
}

TEST(Comparison, Comparison_not_equal)
{
    Sector::Definition sector_def;
    sector_def.location = Sector::ParticleLocation::InsideDetector;
    sector_def.SetMedium(std::make_shared<Medium>(Water()));
    sector_def.SetGeometry(Sphere().create());
    sector_def.scattering_model            = ScatteringFactory::Moliere;
    sector_def.cut_settings                = EnergyCutSettings();
    sector_def.do_continuous_randomization = true;

    std::vector<Sector::Definition> sec_defs;
    sec_defs.push_back(sector_def);

    Propagator prop_a(MuMinusDef::Get(), sec_defs, Sphere().create());
    Propagator prop_b(TauMinusDef::Get(), sec_defs, Sphere().create());

    EXPECT_TRUE(prop_a != prop_b);

    Propagator prop_c(MuMinusDef::Get(), sec_defs, Sphere().create());
    Propagator prop_d(MuMinusDef::Get(), sec_defs, Box().create());

    EXPECT_TRUE(prop_c != prop_d);

    std::vector<Sector::Definition> sec_defs_diff_1;
    sec_defs_diff_1.push_back(sector_def);
    sec_defs_diff_1.push_back(sector_def);

    Propagator prop_e(MuMinusDef::Get(), sec_defs, Sphere().create());
    Propagator prop_f(MuMinusDef::Get(), sec_defs_diff_1, Sphere().create());

    EXPECT_TRUE(prop_e != prop_f);

    // Sector defs differ by the medium
    Sector::Definition sector_def_2;
    sector_def_2.location = Sector::ParticleLocation::InsideDetector;
    sector_def_2.SetMedium(std::make_shared<Medium>(Ice()));
    sector_def_2.SetGeometry(Sphere().create());
    sector_def_2.scattering_model            = ScatteringFactory::Moliere;
    sector_def_2.cut_settings                = EnergyCutSettings();
    sector_def_2.do_continuous_randomization = true;

    std::vector<Sector::Definition> sec_defs_diff_2;
    sec_defs_diff_2.push_back(sector_def_2);

    Propagator prop_g(MuMinusDef::Get(), sec_defs, Sphere().create());
    Propagator prop_h(MuMinusDef::Get(), sec_defs_diff_2, Sphere().create());

    EXPECT_TRUE(prop_g != prop_h);
}

TEST(Assignment, Copyconstructor)
{
    Sector::Definition sector_def;
    sector_def.location = Sector::ParticleLocation::InsideDetector;
    sector_def.SetMedium(std::make_shared<Medium>(Water()));
    sector_def.SetGeometry(Sphere().create());
    sector_def.scattering_model            = ScatteringFactory::Moliere;
    sector_def.cut_settings                = EnergyCutSettings();
    sector_def.do_continuous_randomization = true;

    std::vector<Sector::Definition> sec_defs;
    sec_defs.push_back(sector_def);

    Propagator prop_a(MuMinusDef::Get(), sec_defs, Sphere().create());
    Propagator prop_b = prop_a;

    EXPECT_TRUE(prop_a == prop_b);
}

TEST(Assignment, Copyconstructor2)
{
    Sector::Definition sector_def;
    sector_def.location = Sector::ParticleLocation::InsideDetector;
    sector_def.SetMedium(std::make_shared<Medium>(Water()));
    sector_def.SetGeometry(Sphere().create());
    sector_def.scattering_model            = ScatteringFactory::Moliere;
    sector_def.cut_settings                = EnergyCutSettings();
    sector_def.do_continuous_randomization = true;

    std::vector<Sector::Definition> sec_defs;
    sec_defs.push_back(sector_def);

    Propagator prop_a(MuMinusDef::Get(), sec_defs, Sphere().create());
    Propagator prop_b(prop_a);

    EXPECT_TRUE(prop_a == prop_b);
}

TEST(Propagation, Test_nan)
{
    int statistic = 10;
    int EmaxLog10 = 8;

    // Define Particles
    ParticleDef mu_def = MuMinusDef::Get();
    ParticleDef tau_def = TauMinusDef::Get();
    DynamicData mu(mu_def.particle_type);
    DynamicData tau(tau_def.particle_type);

    // one propagator for each particle definition
    // medium/detector configuration
    Propagator prop_mu(mu_def, "resources/config_ice.json");
    Propagator prop_tau(tau_def, "resources/config_ice.json");

    // Possibility to register propagator in a service
    PropagatorService prop_service;
    prop_service.RegisterPropagator(prop_mu);
    prop_service.RegisterPropagator(prop_tau);

    for (int i = 0; i < statistic; i++)
    {
        // ----------------------------------------------------------------- //
        // Using propagator service
        // ----------------------------------------------------------------- //

        // Set particle properties
        mu.SetEnergy(std::pow(10, EmaxLog10));
        mu.SetPropagatedDistance(0);
        mu.SetPosition(Vector3D(0, 0, 0));
        mu.SetDirection(Vector3D(0, 0, -1));

        tau.SetEnergy(std::pow(10, EmaxLog10));
        tau.SetPropagatedDistance(0);
        tau.SetPosition(Vector3D(0, 0, 0));
        tau.SetDirection(Vector3D(0, 0, -1));

        // Use service to propagate different particle
        Secondaries sec_mu  = prop_service.Propagate(mu_def, mu);
        Secondaries sec_tau = prop_service.Propagate(tau_def, tau);

        // ----------------------------------------------------------------- //
        // Using propagator directly
        // ----------------------------------------------------------------- //

        // Therefor its needed to get the internal created particle first
        // Particle& particle = prop_mu.GetParticle();

        mu.SetEnergy(std::pow(10, EmaxLog10));
        mu.SetPropagatedDistance(0);
        mu.SetPosition(Vector3D(0, 0, 0));
        mu.SetDirection(Vector3D(0, 0, -1));

        Secondaries sec_mu_direct = prop_mu.Propagate(mu);
    }
}

TEST(Propagation, particle_type)
{
    std::string filename = "bin/TestFiles/Propagator_propagation.txt";
	std::ifstream in{filename};
	EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    // Just skip the header
    char firstLine[256];
    in.getline(firstLine, 256);

    int statistic = 10;
    double energy = 1e8;

    in >> statistic >> energy;

    ParticleDef mu_def = MuMinusDef::Get();
    Propagator prop_mu(mu_def, "resources/config_ice.json");
    DynamicData mu(mu_def.particle_type);

    mu.SetEnergy(energy);
    mu.SetPropagatedDistance(0);
    mu.SetPosition(Vector3D(0, 0, 0));
    mu.SetDirection(Vector3D(0, 0, -1));

    std::vector<std::string> names;
    std::vector<double> lengths;
    std::vector<double> energies;
    std::vector<Vector3D> positions;
    std::vector<Vector3D> directions;


    // Buffer to read data in

    std::string name;
    double length     = 0.0;
    double sec_energy = 0.0, x = 0.0, y = 0.0, z = 0.0, dx = 0.0, dy = 0.0, dz = 0.0;

    double error = 1e-2;

    // Read

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);

    for (int i = 0; i < statistic; ++i)
    {
        std::vector<DynamicData> sec_mu_direct = prop_mu.Propagate(mu).GetSecondaries();

        std::string name_new = "";

        for (unsigned int j = 0; j < sec_mu_direct.size(); ++j)
        {
            name_new = sec_mu_direct[j].GetName();

            in >> name >> length >> sec_energy >> x >> y >> z >> dx >> dy >> dz;

            double energy_new  = sec_mu_direct[j].GetEnergy();
            double lenght_new  = sec_mu_direct[j].GetPropagatedDistance();
            Vector3D position  = sec_mu_direct[j].GetPosition();
            Vector3D direction = sec_mu_direct[j].GetDirection();

            EXPECT_TRUE(name == name_new);
            ASSERT_NEAR(length, lenght_new, error * length);
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

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
