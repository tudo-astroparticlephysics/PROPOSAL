
#include "gtest/gtest.h"

#include "PROPOSAL/PROPOSAL.h"

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
    Particle Mu(MuMinusDef::Get());
    Water medium;
    Sphere sphere;
    EnergyCutSettings ecuts;

    Sector::Definition sector_def;
    sector_def.location = Sector::ParticleLocation::InsideDetector;
    sector_def.SetMedium(medium);
    sector_def.SetGeometry(sphere);
    sector_def.scattering_model            = ScatteringFactory::Moliere;
    sector_def.cut_settings                = ecuts;
    sector_def.do_continuous_randomization = true;

    Sector sector(Mu, sector_def);
    Sector sector_2 = Sector(Mu, sector_def);

    EXPECT_TRUE(sector == sector_2);
}

TEST(Comparison, Comparison_not_equal)
{
    Particle Mu(MuMinusDef::Get());
    Particle Tau(TauMinusDef::Get());
    Water medium1;
    Ice medium2;
    Sphere geometry1;
    Cylinder geometry2;
    EnergyCutSettings ecuts1(500, 0.05);
    EnergyCutSettings ecuts2(400, 0.005);

    Sector::Definition sector_def1;
    sector_def1.location = Sector::ParticleLocation::InsideDetector;
    sector_def1.SetMedium(medium1);
    sector_def1.SetGeometry(geometry1);
    sector_def1.scattering_model            = ScatteringFactory::Moliere;
    sector_def1.cut_settings                = ecuts1;
    sector_def1.do_continuous_randomization = true;
    sector_def1.do_exact_time_calculation   = true;
    sector_def1.stopping_decay              = true;

    Sector::Definition sector_def2 = sector_def1;
    sector_def2.location           = Sector::ParticleLocation::InfrontDetector;

    Sector::Definition sector_def3 = sector_def1;
    sector_def3.SetMedium(medium2);

    Sector::Definition sector_def4 = sector_def1;
    sector_def4.SetGeometry(geometry2);

    Sector::Definition sector_def5 = sector_def1;
    sector_def5.scattering_model   = ScatteringFactory::Highland;

    Sector::Definition sector_def6 = sector_def1;
    sector_def6.cut_settings       = ecuts2;

    Sector::Definition sector_def7          = sector_def1;
    sector_def7.do_continuous_randomization = false;

    Sector::Definition sector_def8        = sector_def1;
    sector_def8.do_exact_time_calculation = false;

    Sector::Definition sector_def9 = sector_def1;
    sector_def9.stopping_decay     = false;

    Sector sector(Mu, sector_def1);
    Sector sector_1(Tau, sector_def1);
    Sector sector_2(Mu, sector_def2);
    Sector sector_3(Mu, sector_def3);
    Sector sector_4(Mu, sector_def4);
    Sector sector_5(Mu, sector_def5);
    Sector sector_6(Mu, sector_def6);
    Sector sector_7(Mu, sector_def7);
    Sector sector_8(Mu, sector_def8);
    Sector sector_9(Mu, sector_def9);

    EXPECT_TRUE(sector != sector_1);
    EXPECT_TRUE(sector != sector_2);
    EXPECT_TRUE(sector != sector_3);
    EXPECT_TRUE(sector != sector_4);
    EXPECT_TRUE(sector != sector_5);
    EXPECT_TRUE(sector != sector_6);
    EXPECT_TRUE(sector != sector_7);
    EXPECT_TRUE(sector != sector_8);
    EXPECT_TRUE(sector != sector_9);
}

TEST(Assignment, Copyconstructor)
{
    Particle mu = Particle(MuMinusDef::Get());
    Water water(1.0);
    Sphere geometry(Vector3D(0, 0, 0), 1000, 0);
    EnergyCutSettings ecuts(500, 0.05);

    Sector::Definition sector_def;
    sector_def.location = Sector::ParticleLocation::InsideDetector;
    sector_def.SetMedium(water);
    sector_def.SetGeometry(geometry);
    sector_def.scattering_model = ScatteringFactory::Moliere;
    sector_def.cut_settings     = ecuts;

    Sector sector_1(mu, sector_def);
    Sector sector_2 = sector_1;
    EXPECT_TRUE(sector_1 == sector_2);
}

TEST(Assignment, Copyconstructor2)
{
    Particle mu = Particle(MuMinusDef::Get());
    Water water(1.0);
    Sphere geometry(Vector3D(), 1000, 0);
    EnergyCutSettings ecuts(500, 0.05);

    Sector::Definition sector_def;
    sector_def.location = Sector::ParticleLocation::InsideDetector;
    sector_def.SetMedium(water);
    sector_def.SetGeometry(geometry);
    sector_def.scattering_model = ScatteringFactory::Moliere;
    sector_def.cut_settings     = ecuts;

    Sector sector_1(mu, sector_def);
    Sector sector_2(sector_1);
    EXPECT_TRUE(sector_1 == sector_2);
}

TEST(Sector, Propagate)
{
    std::ifstream in;
    std::string filename = "bin/TestFiles/Sector_propagate.txt";
    in.open(filename.c_str());

    if (!in.good())
    {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    std::string particleName;
    std::string mediumName;

    Utility::Definition utility_def;

    double energyTillStochastic_calc;
    double energyTillStochastic_stored;
    double stochasticLoss_calc;
    double stochasticLoss_stored;
    double energy_final_calc;
    double energy_final_stored;
    double energy_previous = -1;
    double energy_init;
    double distance;
    double ecut, vcut;

    bool first_line = true;

    std::cout.precision(16);

    RandomGenerator::Get().SetSeed(1234);

    while (in.good())
    {
        if (first_line)
        {
            in >> particleName >> mediumName >> ecut >> vcut >> energy_init >> energyTillStochastic_stored >> stochasticLoss_stored >> energy_final_stored >> distance;
            first_line = false;
        }

        energy_previous   = -1;
        Particle particle = Particle(getParticleDef(particleName));
        particle.SetEnergy(energy_init);
        particle.SetDirection(Vector3D(1, 0, 0));
        Medium* medium = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        Sector::Definition sector_def;
        sector_def.SetMedium(*medium);
        sector_def.cut_settings = ecuts;

        Sector sector(particle, sector_def, InterpolationDef());

        while (energy_previous < energy_init)
        {
            energy_previous = energy_init;
            particle.SetEnergy(energy_init);
            particle.SetDirection(Vector3D(1, 0, 0));

            energyTillStochastic_calc = sector.CalculateEnergyTillStochastic(energy_init).first;
            stochasticLoss_calc = sector.MakeStochasticLoss(energy_init).first;
            energy_final_calc = sector.Propagate(distance);

            ASSERT_NEAR(energyTillStochastic_calc, energyTillStochastic_stored, std::abs(1e-3 * energyTillStochastic_calc));
            ASSERT_NEAR(stochasticLoss_calc, stochasticLoss_stored, std::abs(1e-3 * stochasticLoss_calc));
            ASSERT_NEAR(energy_final_calc, energy_final_stored, std::abs(1e-3 * energy_final_calc));

            in >> particleName >> mediumName >> ecut >> vcut >> energy_init >> energyTillStochastic_stored >> stochasticLoss_stored >> energy_final_stored >> distance;
        }

        delete medium;
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
