
#include "gtest/gtest.h"

#include "PROPOSAL/PROPOSAL.h"
#include <string>
using namespace PROPOSAL;

/* std::string PATH_TO_TABLES = "~/.local/share/PROPOSAL/tables"; */
std::string PATH_TO_TABLES = "";

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
    Particle Mu(MuMinusDef::Get());
    Water medium;
    Sphere sphere;
    EnergyCutSettings ecuts;

    Sector::Definition sector_def;
    sector_def.location = Sector::ParticleLocation::InsideDetector;
    sector_def.SetMedium(medium);
    sector_def.SetGeometry(sphere);
    sector_def.scattering_model = ScatteringFactory::Moliere;
    sector_def.cut_settings = ecuts;
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
    sector_def.cut_settings = ecuts;

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
    sector_def.cut_settings = ecuts;

    Sector sector_1(mu, sector_def);
    Sector sector_2(sector_1);
    EXPECT_TRUE(sector_1 == sector_2);
}

TEST(Sector, Continous)
{
    std::ifstream in;
    std::string filename = "bin/TestFiles/Sector_ContinousLoss.txt";
    in.open(filename.c_str());

    if (!in.good()) {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

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
    Particle* particle = new Particle(MuMinusDef::Get());
    Medium* medium = new Medium(sector_def.GetMedium());
    EnergyCutSettings* cuts = new EnergyCutSettings(sector_def.cut_settings);
    Sector* sector = new Sector(*particle, sector_def, inter_def);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        if (particle->GetParticleDef().name != particleName
            || cuts->GetEcut() != ecut || cuts->GetVcut() != vcut
            || medium->GetName() != mediumName) {
            delete particle, sector, medium, cuts;

            cuts = new EnergyCutSettings(ecut, vcut);
            medium = MediumFactory::Get().CreateMedium(mediumName);
            particle = new Particle(getParticleDef(particleName));

            sector_def.SetMedium(*medium);
            sector_def.cut_settings = *cuts;

            sector = new Sector(*particle, sector_def, inter_def);
        }

        while (ss >> energy) {
            std::pair<double, double> energies
                = sector->CalculateEnergyTillStochastic(initial_energy);
            double energy_calc = std::max(energies.first, energies.second);
            ASSERT_NEAR(energy_calc, energy, std::abs(1e-3 * energy_calc));
            initial_energy = energy;
        }
    }

}

TEST(Sector, Stochastic)
{
    std::ifstream in;
    std::string filename = "bin/TestFiles/Sector_StochasticLoss.txt";
    in.open(filename.c_str());

    if (!in.good()) {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

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
    Particle* particle = new Particle(MuMinusDef::Get());
    Medium* medium = new Medium(sector_def.GetMedium());
    EnergyCutSettings* cuts = new EnergyCutSettings(sector_def.cut_settings);
    Sector* sector = new Sector(*particle, sector_def, inter_def);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        if (particle->GetParticleDef().name != particleName
            || cuts->GetEcut() != ecut || cuts->GetVcut() != vcut
            || medium->GetName() != mediumName) {
            delete particle, sector, medium, cuts;

            cuts = new EnergyCutSettings(ecut, vcut);
            medium = MediumFactory::Get().CreateMedium(mediumName);
            particle = new Particle(getParticleDef(particleName));

            sector_def.SetMedium(*medium);
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
    std::ifstream in;
    std::string filename = "bin/TestFiles/Sector_Energy_Distance.txt";
    in.open(filename.c_str());

    if (!in.good()) {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

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
    Particle* particle = new Particle(MuMinusDef::Get());
    Medium* medium = new Medium(sector_def.GetMedium());
    EnergyCutSettings* cuts = new EnergyCutSettings(sector_def.cut_settings);
    Sector* sector = new Sector(*particle, sector_def, inter_def);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        if (particle->GetParticleDef().name != particleName
            || cuts->GetEcut() != ecut || cuts->GetVcut() != vcut
            || medium->GetName() != mediumName) {
            delete particle, sector, medium, cuts;

            cuts = new EnergyCutSettings(ecut, vcut);
            medium = MediumFactory::Get().CreateMedium(mediumName);
            particle = new Particle(getParticleDef(particleName));

            sector_def.SetMedium(*medium);
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
    std::ifstream in;
    std::string filename = "bin/TestFiles/Sector_Propagate.txt";
    in.open(filename.c_str());

    if (!in.good()) {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

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
    Particle* particle = new Particle(MuMinusDef::Get());
    Medium* medium = new Medium(sector_def.GetMedium());
    EnergyCutSettings* cuts = new EnergyCutSettings(sector_def.cut_settings);
    Sector* sector = new Sector(*particle, sector_def, inter_def);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        if (particle->GetParticleDef().name != particleName
            || cuts->GetEcut() != ecut || cuts->GetVcut() != vcut
            || medium->GetName() != mediumName) {
            delete particle, sector, medium, cuts;

            cuts = new EnergyCutSettings(ecut, vcut);
            medium = MediumFactory::Get().CreateMedium(mediumName);
            particle = new Particle(getParticleDef(particleName));

            sector_def.SetMedium(*medium);
            sector_def.cut_settings = *cuts;

            sector = new Sector(*particle, sector_def, inter_def);
        }

        while (ss >> produced_particle >> energy) {
            particle->SetDirection(Vector3D(0, 0, -1));
            particle->SetPosition(Vector3D(0, 0, 0));
            particle->SetEnergy(initial_energy);

            std::pair<double, Secondaries> aux = sector->Propagate(1000);
            int produced_particle = aux.second.GetNumberOfParticles();
            ASSERT_NEAR(aux.first, energy, std::abs(1e-3 * aux.first));
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
