
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/Factories/WeakInteractionFactory.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSALTestUtilities/TestFilesHandling.h"
#include "PROPOSAL/secondaries/parametrization/weakinteraction/WeakCooperSarkarMertsch.h"

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

TEST(WeakInteraction, Test_of_dNdx)
{
    std::ifstream in;
    getTestFile("Weak_dNdx.txt", in);

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

        auto cross = make_weakinteraction(particle_def, *medium, false, config);

        dNdx_new = cross->CalculatedNdx(energy);

        EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);
    }
}

TEST(WeakInteraction, Test_Stochastic_Loss)
{
    std::ifstream in;
    getTestFile("Weak_e.txt", in);

    std::string particleName;
    std::string mediumName;
    double multiplier;
    double energy;
    std::string parametrization;
    double rnd1;
    double rnd2;
    double v_stored;

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> multiplier >> energy
            >> parametrization >> rnd1 >> rnd2 >> v_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto secondaries = secondaries::WeakCooperSarkarMertsch(particle_def, *medium);

        auto components = medium->GetComponents();
        auto comp = components.at(int(rnd2 * components.size()));

        auto v = secondaries.CalculateRelativeLoss(energy, rnd1, comp);
        EXPECT_NEAR(v, v_stored, v_stored * 1e-3);
    }
}

TEST(WeakInteraction, Test_of_dNdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Weak_dNdx.txt", in);

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

        auto cross = make_weakinteraction(particle_def, *medium, true, config);

        dNdx_new = cross->CalculatedNdx(energy);

        EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
    }
}
int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
