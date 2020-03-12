
#include "gtest/gtest.h"
#include <algorithm>

#include "PROPOSAL/PROPOSAL.h"
#include <string>
using namespace PROPOSAL;

std::string PATH_TO_TABLES = "~/.local/share/PROPOSAL/tables";
/* std::string PATH_TO_TABLES = ""; */

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "MuMinus") {
        return MuMinusDef::Get();
    } else if (name == "TauMinus") {
        return TauMinusDef::Get();
    } else {
        return EMinusDef::Get();
    }
}

TEST(Comparison, Comparison_equal)
{
    ParticleDef mu_def = MuMinusDef::Get();
    auto medium  = std::make_shared<Water>();
    auto sphere = Sphere().create();
    EnergyCutSettings ecuts;

    Sector::Definition sector_def;
    sector_def.location = Sector::ParticleLocation::InsideDetector;
    sector_def.SetMedium(medium);
    sector_def.SetGeometry(sphere);
    sector_def.scattering_model = ScatteringFactory::Moliere;
    sector_def.cut_settings = ecuts;
    sector_def.do_continuous_randomization = true;

    Sector sector(mu_def, sector_def);
    Sector sector_2 = Sector(mu_def, sector_def);

    EXPECT_TRUE(sector == sector_2);
}

TEST(Comparison, Comparison_not_equal)
{
    ParticleDef mu_def = MuMinusDef::Get();
    ParticleDef tau_def = TauMinusDef::Get();
    auto medium1 = std::make_shared<Water>();
    auto medium2 = std::make_shared<Ice>();
    auto geometry1 = Sphere().create();
    auto geometry2 = Cylinder().create();
    EnergyCutSettings ecuts1(500, 0.05);
    EnergyCutSettings ecuts2(400, 0.005);

    Sector::Definition sector_def1;
    sector_def1.location = Sector::ParticleLocation::InsideDetector;
    sector_def1.SetMedium(medium1);
    sector_def1.SetGeometry(geometry1);
    sector_def1.scattering_model = ScatteringFactory::Moliere;
    sector_def1.cut_settings = ecuts1;
    sector_def1.do_continuous_randomization = true;
    sector_def1.do_exact_time_calculation = true;
    sector_def1.stopping_decay = true;

    Sector::Definition sector_def2 = sector_def1;
    sector_def2.location = Sector::ParticleLocation::InfrontDetector;

    Sector::Definition sector_def3 = sector_def1;
    sector_def3.SetMedium(medium2);

    Sector::Definition sector_def4 = sector_def1;
    sector_def4.SetGeometry(geometry2);

    Sector::Definition sector_def5 = sector_def1;
    sector_def5.scattering_model = ScatteringFactory::Highland;

    Sector::Definition sector_def6 = sector_def1;
    sector_def6.cut_settings = ecuts2;

    Sector::Definition sector_def7 = sector_def1;
    sector_def7.do_continuous_randomization = false;

    Sector::Definition sector_def8 = sector_def1;
    sector_def8.do_exact_time_calculation = false;

    Sector::Definition sector_def9 = sector_def1;
    sector_def9.stopping_decay = false;

    Sector sector(mu_def, sector_def1);
    Sector sector_1(tau_def, sector_def1);
    Sector sector_2(mu_def, sector_def2);
    Sector sector_3(mu_def, sector_def3);
    Sector sector_4(mu_def, sector_def4);
    Sector sector_5(mu_def, sector_def5);
    Sector sector_6(mu_def, sector_def6);
    Sector sector_7(mu_def, sector_def7);
    Sector sector_8(mu_def, sector_def8);
    Sector sector_9(mu_def, sector_def9);

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
    ParticleDef mu_def = MuMinusDef::Get();
    auto water = std::make_shared<Water>();
    auto geometry = Sphere(Vector3D(0, 0, 0), 1000, 0).create();
    EnergyCutSettings ecuts(500, 0.05);

    Sector::Definition sector_def;
    sector_def.location = Sector::ParticleLocation::InsideDetector;
    sector_def.SetMedium(water);
    sector_def.SetGeometry(geometry);
    sector_def.scattering_model = ScatteringFactory::Moliere;
    sector_def.cut_settings = ecuts;

    Sector sector_1(mu_def, sector_def);
    Sector sector_2 = sector_1;
    EXPECT_TRUE(sector_1 == sector_2);
}

TEST(Assignment, Copyconstructor2)
{
    ParticleDef mu = MuMinusDef::Get();
    auto water = std::make_shared<Water>();
    auto geometry = Sphere(Vector3D(), 1000, 0).create();
    EnergyCutSettings ecuts(500, 0.05);

    Sector::Definition sector_def;
    sector_def.location = Sector::ParticleLocation::InsideDetector;
    sector_def.SetMedium(water);
    sector_def.SetGeometry(geometry);
    sector_def.scattering_model = ScatteringFactory::Moliere;
    sector_def.cut_settings = ecuts;

    Sector sector_1(mu, sector_def);
    Sector sector_2(sector_1);
    EXPECT_TRUE(sector_1 == sector_2);
}

TEST(Sector, Continuous)
{
    std::string filename = "bin/TestFiles/Sector_ContinousLoss.txt";
	std::ifstream in{filename};
	EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);
    std::string particleName;
    std::string mediumName;
    double ecut, vcut;
    double energy, initial_energy;
    InterpolationDef inter_def;
    inter_def.path_to_tables = PATH_TO_TABLES;
    inter_def.path_to_tables_readonly = PATH_TO_TABLES;

    Sector::Definition sector_def;
    ParticleDef* particle = new ParticleDef(MuMinusDef::Get());
    std::shared_ptr<const Medium> medium = sector_def.GetMedium();
    EnergyCutSettings* cuts = new EnergyCutSettings(sector_def.cut_settings);
    Sector* sector = new Sector(*particle, sector_def, inter_def);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        if (particle->name != particleName || cuts->GetEcut() != ecut
            || cuts->GetVcut() != vcut || medium->GetName() != mediumName) {
            delete particle, sector, medium, cuts;

            cuts = new EnergyCutSettings(ecut, vcut);
            std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);
            particle = new ParticleDef(getParticleDef(particleName));

            sector_def.SetMedium(medium);
            sector_def.cut_settings = *cuts;

            sector = new Sector(*particle, sector_def, inter_def);
        }

        while (ss >> energy) {
            double rndd = RandomGenerator::Get().RandomDouble();
            double decay_energy = sector->EnergyDecay(initial_energy, rndd);
            double rndi = RandomGenerator::Get().RandomDouble();
            double inter_energy
                = sector->EnergyInteraction(initial_energy, rndi);
            double energy_calc = std::max(decay_energy, inter_energy);
            ASSERT_NEAR(energy_calc, energy, std::abs(1e-3 * energy_calc));
            initial_energy = energy;
        }
    }
}

TEST(Sector, Stochastic)
{
    std::string filename = "bin/TestFiles/Sector_StochasticLoss.txt";
	std::ifstream in{filename};
	EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);
    std::string particleName;
    std::string mediumName;
    double ecut, vcut;
    double energy, initial_energy, rnd;
    int interaction_type;
    InterpolationDef inter_def;
    inter_def.path_to_tables = PATH_TO_TABLES;
    inter_def.path_to_tables_readonly = PATH_TO_TABLES;

    Sector::Definition sector_def;
    ParticleDef* particle = new ParticleDef(MuMinusDef::Get());
    std::shared_ptr<const Medium> medium=sector_def.GetMedium();
    EnergyCutSettings* cuts = new EnergyCutSettings(sector_def.cut_settings);
    Sector* sector = new Sector(*particle, sector_def, inter_def);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        if (particle->name != particleName || cuts->GetEcut() != ecut
            || cuts->GetVcut() != vcut || medium->GetName() != mediumName) {
            delete particle, sector, medium, cuts;

            cuts = new EnergyCutSettings(ecut, vcut);
            std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);
            particle = new ParticleDef(getParticleDef(particleName));

            sector_def.SetMedium(medium);
            sector_def.cut_settings = *cuts;

            sector = new Sector(*particle, sector_def, inter_def);
        }

        while (ss >> energy >> interaction_type >> rnd) {
            std::pair<double, int> loss =sector->MakeStochasticLoss(initial_energy);
            double energy_calc = initial_energy - loss.first;
            double random = RandomGenerator::Get().RandomDouble();
            ASSERT_NEAR(random, rnd, std::abs(1e-3 * energy_calc));
            ASSERT_NEAR(loss.second, interaction_type, std::abs(1e-3 * energy_calc));
            ASSERT_NEAR(energy_calc, energy, std::abs(1e-3 * energy_calc));
            initial_energy = energy;
        }
    }
}

TEST(Sector, EnergyDisplacement)
{
    std::string filename = "bin/TestFiles/Sector_Energy_Distance.txt";
	std::ifstream in{filename};
	EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);
    std::string particleName;
    std::string mediumName;
    double ecut, vcut;
    double displacement, energy, initial_energy;
    InterpolationDef inter_def;
    inter_def.path_to_tables = PATH_TO_TABLES;
    inter_def.path_to_tables_readonly = PATH_TO_TABLES;

    Sector::Definition sector_def;
    ParticleDef* particle = new ParticleDef(MuMinusDef::Get());
    std::shared_ptr<const Medium> medium(sector_def.GetMedium());
    EnergyCutSettings* cuts = new EnergyCutSettings(sector_def.cut_settings);
    Sector* sector = new Sector(*particle, sector_def, inter_def);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        if (particle->name != particleName || cuts->GetEcut() != ecut
            || cuts->GetVcut() != vcut || medium->GetName() != mediumName) {
            delete particle, sector, medium, cuts;

            cuts = new EnergyCutSettings(ecut, vcut);
            std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);
            particle = new ParticleDef(getParticleDef(particleName));

            sector_def.SetMedium(medium);
            sector_def.cut_settings = *cuts;

            sector = new Sector(*particle, sector_def, inter_def);
        }

        while (ss >> displacement >> energy) {
            double energy_calc
                = sector->EnergyDistance(initial_energy, displacement);
            ASSERT_NEAR(energy_calc, energy, std::abs(1e-3 * energy_calc));
        }
    }
}

TEST(Sector, Propagate)
{
    std::string filename = "bin/TestFiles/Sector_Propagate.txt";
	std::ifstream in{filename};
	EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);
    std::string particleName;
    std::string mediumName;
    double ecut, vcut;
    double energy, initial_energy;
    int produced_particle;
    InterpolationDef inter_def;
    inter_def.path_to_tables = PATH_TO_TABLES;
    inter_def.path_to_tables_readonly = PATH_TO_TABLES;

    Sector::Definition sector_def;
    ParticleDef* particle = new ParticleDef(MuMinusDef::Get());
    std::shared_ptr<const Medium> medium(sector_def.GetMedium());
    EnergyCutSettings* cuts = new EnergyCutSettings(sector_def.cut_settings);
    Sector* sector = new Sector(*particle, sector_def, inter_def);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        if (particle->name != particleName || cuts->GetEcut() != ecut
            || cuts->GetVcut() != vcut || medium->GetName() != mediumName) {
            delete particle, sector, medium, cuts;

            cuts = new EnergyCutSettings(ecut, vcut);
            std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);
            particle = new ParticleDef(getParticleDef(particleName));

            sector_def.SetMedium(medium);
            sector_def.cut_settings = *cuts;

            sector = new Sector(*particle, sector_def, inter_def);
        }

        while (ss >> produced_particle >> energy) {
            DynamicData p_condition;
            p_condition.SetDirection(Vector3D(0, 0, -1));
            p_condition.SetPosition(Vector3D(0, 0, 0));
            p_condition.SetEnergy(initial_energy);

            Secondaries secondaries = sector->Propagate(p_condition, 1000, 0);
            int produced_particle = secondaries.GetNumberOfParticles();

            double energy_calc = secondaries.GetSecondaries().back().GetEnergy();
            DynamicData last_condition = secondaries.GetSecondaries().back();
            if (last_condition.GetType() != static_cast<int>(InteractionType::Decay)) {
                ASSERT_NEAR(last_condition.GetEnergy(), energy, std::abs(1e-3 * energy_calc));
            } else {
                ASSERT_NEAR(-last_condition.GetPropagatedDistance(), energy, std::abs(1e-3 * energy_calc));
            }
            /* double energy_calc = secondaries.GetSecondaries().back().GetEnergy(); */
            /* ASSERT_NEAR(energy_calc, energy, std::abs(1e-3 * energy_calc)); */
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
