#include "gtest/gtest.h"
#include <nlohmann/json.hpp>

#include "PROPOSAL/crosssection/Factories/PhotoMuPairProductionFactory.h"
#include "PROPOSAL/secondaries/parametrization/photomupairproduction/PhotoMuPairProductionBurkhardtKelnerKokoulin.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSALTestUtilities/TestFilesHandling.h"

using namespace PROPOSAL;

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "Gamma") {
        return GammaDef();
    } else {
        EXPECT_TRUE(false);
        return MuMinusDef();
    }
}

TEST(PhotoMuPair, Test_of_dNdx)
{
    std::ifstream in;
    getTestFile("PhotoMuPair_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    double multiplier;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> multiplier >> energy
        >> parametrization >> dNdx_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_photomupairproduction(particle_def, *medium, false, config);

        dNdx_new = cross->CalculatedNdx(energy);
        EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);

        auto dEdx = cross->CalculatedEdx(energy);
        EXPECT_EQ(dEdx, 0); // only-stochastic interaction

        auto dE2dx = cross->CalculatedE2dx(energy);
        EXPECT_EQ(dE2dx, 0); // only-stochastic interaction

        auto v_loss = cross->CalculateStochasticLoss(0, energy, dNdx_new);
        EXPECT_EQ(v_loss, 1); // fatal interaction
    }
}

TEST(PhotoMuPair, Test_of_dNdx_Interpolant)
{
    std::ifstream in;
    getTestFile("PhotoMuPair_dNdx.txt", in);

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

        auto cross = make_photomupairproduction(particle_def, *medium, true, config);

        dNdx_new = cross->CalculatedNdx(energy);
        EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);

        auto dEdx = cross->CalculatedEdx(energy);
        EXPECT_EQ(dEdx, 0); // only-stochastic interaction

        auto dE2dx = cross->CalculatedE2dx(energy);
        EXPECT_EQ(dE2dx, 0); // only-stochastic interaction

        auto v_loss = cross->CalculateStochasticLoss(0, energy, dNdx_new);
        EXPECT_EQ(v_loss, 1); // fatal interaction

    }
}

TEST(PhotoMuPair, CalculateX)
{
    std::ifstream in;
    getTestFile("PhotoMuPair_x.txt", in);

    std::string particleName;
    std::string mediumName;
    double rnd1;
    double rnd2;
    double energy;
    double x_stored;
    double x;
    std::string parametrization;

    while (in >> particleName >> mediumName >> energy >> parametrization >> x_stored
        >> rnd1 >> rnd2) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        std::unique_ptr<secondaries::PhotoMuPairProduction> param;
        if (parametrization == "BurkhardtKelnerKokoulin")
            param = std::make_unique<secondaries::PhotoMuPairProductionBurkhardtKelnerKokoulin>(particle_def, *medium);
        else
            throw std::logic_error("Unknown parametrization in PhotoMuPair_TEST");

        auto components = medium->GetComponents();
        auto comp = components.at(int(rnd2 * components.size()));

        x = param->Calculatex(energy, rnd1, comp);
        EXPECT_NEAR(x, x_stored, 1e-10 * x_stored);
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
