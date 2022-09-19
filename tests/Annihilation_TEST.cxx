
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/Factories/AnnihilationFactory.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSALTestUtilities/TestFilesHandling.h"

using namespace PROPOSAL;

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "MuMinus") {
        return MuMinusDef();
    } else if (name == "TauMinus") {
        return TauMinusDef();
    } else if (name == "EMinus") {
        return EMinusDef();
    } else if (name == "MuPlus") {
        return MuPlusDef();
    } else if (name == "TauPlus") {
        return TauPlusDef();
    } else if (name == "EPlus") {
        return EPlusDef();
    } else {
        return MuMinusDef();
    }
}

TEST(Annihilation, Test_of_dEdx)
{
    ParticleDef particle_def = EPlusDef();
    auto medium = Air();
    nlohmann::json config;
    config["parametrization"] = "Heitler";
    auto cross_integral = make_annihilation(particle_def, medium, false, config);
    auto cross_interpol = make_annihilation(particle_def, medium, false, config);

    EXPECT_EQ(cross_integral->CalculatedEdx(1e3), 0);
    EXPECT_EQ(cross_integral->CalculatedE2dx(1e3), 0);

    EXPECT_EQ(cross_interpol->CalculatedEdx(1e3), 0);
    EXPECT_EQ(cross_interpol->CalculatedE2dx(1e3), 0);
}

TEST(Annihilation, Test_of_dNdx)
{
    std::ifstream in;
    getTestFile("Anni_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    double multiplier;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> multiplier >> energy
        >> parametrization >> dNdx_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_annihilation(particle_def, *medium, false, config);

        dNdx_new = cross->CalculatedNdx(energy);
        EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-4 * dNdx_stored);
    }
}

TEST(Annihilation, Test_Stochastic_Loss)
{
    std::ifstream in;
    getTestFile("Anni_e.txt", in);

    std::string particleName;
    std::string mediumName;
    double multiplier;
    double energy;
    std::string parametrization;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> multiplier >> energy
        >> parametrization >> rnd1 >> rnd2 >> stochastic_loss_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_annihilation(particle_def, *medium, false, config);

        auto components = medium->GetComponents();
        auto comp = components.at(int(rnd2 * components.size()));

        double dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());

        auto stochastic_loss = cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp);
        EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-6 * stochastic_loss_stored);
        EXPECT_EQ(stochastic_loss, 1.); // we always loose all energy
    }
}

TEST(Annihilation, Test_of_dNdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Anni_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    double multiplier;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> multiplier >> energy
        >> parametrization >> dNdx_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_annihilation(particle_def, *medium, true, config);

        dNdx_new = cross->CalculatedNdx(energy);
        EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-4 * dNdx_stored);
    }
}

TEST(Annihilation, Test_of_e_Interpolant)
{
    std::ifstream in;
    getTestFile("Anni_e.txt", in);

    std::string particleName;
    std::string mediumName;
    double multiplier;
    double energy;
    std::string parametrization;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;

    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> multiplier >> energy
        >> parametrization >> rnd1 >> rnd2 >> stochastic_loss_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_annihilation(particle_def, *medium, false, config);

        auto components = medium->GetComponents();
        auto comp = components.at(int(rnd2 * components.size()));

        double dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());

        auto stochastic_loss = cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp);
        EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-6 * stochastic_loss_stored);
        EXPECT_EQ(stochastic_loss, 1.); // we always loose all energy
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
