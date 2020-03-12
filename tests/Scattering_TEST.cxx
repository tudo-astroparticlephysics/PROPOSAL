
#include <fstream>

#include "gtest/gtest.h"

#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/scattering/ScatteringFactory.h"
#include "PROPOSAL/scattering/ScatteringHighland.h"
#include "PROPOSAL/scattering/ScatteringHighlandIntegral.h"
#include "PROPOSAL/scattering/ScatteringMoliere.h"
#include "PROPOSAL/scattering/ScatteringNoScattering.h"

#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"

using namespace PROPOSAL;

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "MuMinus")
    {
        return MuMinusDef::Get();
    } else if (name == "TauMinus")
    {
        return TauMinusDef::Get();
    } else
    {
        return EMinusDef::Get();
    }
}

TEST(Comparison, Comparison_equal)
{
    ParticleDef mu = MuMinusDef::Get();
    auto water = std::make_shared<const Water>(1.0);

    Scattering* noScat1 = new ScatteringNoScattering(mu, water);
    ScatteringNoScattering noScat2(mu, water);

    EXPECT_TRUE(*noScat1 == noScat2);

    Scattering* moliere1 = new ScatteringMoliere(mu, water);
    ScatteringMoliere moliere2(mu, water);

    EXPECT_TRUE(*moliere1 == moliere2);

    Scattering* high1 = new ScatteringHighland(mu, water);
    ScatteringHighland high2(mu, water);

    EXPECT_TRUE(*high1 == high2);

    EnergyCutSettings ecuts;
    Utility::Definition utility_defs;
    Utility utils(MuMinusDef::Get(), water, ecuts, utility_defs);

    Scattering* highInt1 = new ScatteringHighlandIntegral(mu, utils);
    ScatteringHighlandIntegral highInt2(mu, utils);

    EXPECT_TRUE(*highInt1 == highInt2);
}

TEST(Comparison, Comparison_not_equal)
{
    ParticleDef mu  = MuMinusDef::Get();
    ParticleDef tau = TauMinusDef::Get();
    auto water = std::make_shared<const Water>(1.0);
    auto ice = std::make_shared<const Ice>();

    ScatteringNoScattering noScat1(mu, water);
    ScatteringNoScattering noScat2(tau, water);
    ScatteringNoScattering noScat3(mu, ice);

    EXPECT_TRUE(noScat1 != noScat2);
    EXPECT_TRUE(noScat1 != noScat3);

    ScatteringMoliere moliere1(mu, water);
    ScatteringMoliere moliere2(tau, water);
    ScatteringMoliere moliere3(mu, ice);

    EXPECT_TRUE(moliere1 != moliere2);
    EXPECT_TRUE(moliere1 != moliere3);

    ScatteringHighland high1(mu, water);
    ScatteringHighland high2(tau, water);
    ScatteringHighland high3(mu, ice);

    EXPECT_TRUE(high1 != high2);
    EXPECT_TRUE(high1 != high3);

    EnergyCutSettings ecuts1;
    EnergyCutSettings ecuts2(200, 0.01);
    Utility::Definition utility_defs;
    Utility utils1(MuMinusDef::Get(), water, ecuts1, utility_defs);
    Utility utils2(TauMinusDef::Get(), water, ecuts1, utility_defs);
    Utility utils3(MuMinusDef::Get(), ice, ecuts1, utility_defs);
    Utility utils4(MuMinusDef::Get(), water, ecuts2, utility_defs);

    ScatteringHighlandIntegral highInt1(mu, utils1);
    ScatteringHighlandIntegral highInt2(tau, utils2);
    ScatteringHighlandIntegral highInt3(mu, utils3);
    ScatteringHighlandIntegral highInt4(mu, utils4);

    EXPECT_TRUE(highInt1 != highInt2);
    EXPECT_TRUE(highInt1 != highInt3);
    EXPECT_TRUE(highInt1 != highInt4);
}

TEST(Assignment, Copyconstructor)
{
    ParticleDef mu = MuMinusDef::Get();
    auto water = std::make_shared<const Water>(1.0);

    ScatteringMoliere moliere1(mu, water);
    ScatteringMoliere moliere2 = moliere1;
    EXPECT_TRUE(moliere1 == moliere2);
}

TEST(Assignment, Copyconstructor2)
{
    ParticleDef mu = MuMinusDef::Get();
    auto water = std::make_shared<const Water>(1.0);

    ScatteringMoliere moliere1(mu, water);
    ScatteringMoliere moliere2(moliere1);
    EXPECT_TRUE(moliere1 == moliere2);
}

TEST(Scattering, Scatter)
{
    std::string filename = "bin/TestFiles/Scattering_scatter.txt";
	std::ifstream in{filename};
	EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::string particleName;
    std::string mediumName;
    std::string parametrization;

    double energy_init, energy_final, distance;
    double rnd1, rnd2, rnd3, rnd4;
    double energy_previous = -1;
    double ecut, vcut;
    Vector3D position_init  = Vector3D(0, 0, 0);
    Vector3D direction_init = Vector3D(1, 0, 0);
    Vector3D position_out;
    Vector3D direction_out;
    direction_init.CalculateSphericalCoordinates();
    double x_f, y_f, z_f;
    double radius_f, phi_f, theta_f;

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);
    double error    = 1e-3;
    bool first_line = true;

    while (in.good())
    {
        if (first_line)
        {
            in >> particleName >> mediumName >> parametrization >> ecut >> vcut >> energy_init >> energy_final >>
                distance >> rnd1 >> rnd2 >> rnd3 >> rnd4 >> x_f >> y_f >> z_f >> radius_f >> phi_f >> theta_f;

            first_line = false;
        }

        energy_previous = -1;

        ParticleDef particle_def = getParticleDef(particleName);

        std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        Scattering* scattering = NULL;

        if (parametrization == "HighlandIntegral")
        {
            Utility utility(particle_def, medium, ecuts, Utility::Definition(), InterpolationDef());
            scattering = ScatteringFactory::Get().CreateScattering(parametrization, particle_def, utility, InterpolationDef());
        }
        else
        {
            Utility utility(particle_def, medium, ecuts, Utility::Definition());
            scattering = ScatteringFactory::Get().CreateScattering(parametrization, particle_def, utility, InterpolationDef());
        }

        while (energy_previous < energy_init)
        {
            energy_previous = energy_init;

            Directions directions = scattering->Scatter(distance,
                                                        energy_init,
                                                        energy_final,
                                                        position_init,
                                                        direction_init,
                                                        rnd1, rnd2, rnd3, rnd4);
            position_out = position_init + distance * directions.u_;
            direction_out = directions.n_i_;

            ASSERT_NEAR(position_out.GetX(), x_f, std::abs(error * x_f));
            ASSERT_NEAR(position_out.GetY(), y_f, std::abs(error * y_f));
            ASSERT_NEAR(position_out.GetZ(), z_f, std::abs(error * z_f));

            ASSERT_NEAR(direction_out.GetRadius(), radius_f, std::abs(error * radius_f));
            ASSERT_NEAR(direction_out.GetPhi(), phi_f, std::abs(error * phi_f));
            ASSERT_NEAR(direction_out.GetTheta(), theta_f, std::abs(error * theta_f));

            in >> particleName >> mediumName >> parametrization >> ecut >> vcut >> energy_init >> energy_final >>
                distance >> rnd1 >> rnd2 >> rnd3 >> rnd4 >> x_f >> y_f >> z_f >> radius_f >> phi_f >> theta_f;
        }

        delete scattering;
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
