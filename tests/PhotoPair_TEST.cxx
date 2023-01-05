
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/Factories/PhotoPairProductionFactory.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
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

TEST(PhotoPair, Test_of_dNdx)
{
    std::ifstream in;
    getTestFile("PhotoPair_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    double multiplier;
    bool lpm;
    double density_correction;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> multiplier >> lpm >> density_correction >> energy
        >> parametrization >> dNdx_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["lpm"] = lpm;

        auto cross
            = make_photopairproduction(particle_def, *medium, false, config, density_correction);

        dNdx_new = cross->CalculatedNdx(energy);
        EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-5 * dNdx_stored);
    }
}

TEST(PhotoPair, Test_of_dNdx_Interpolant)
{
    std::ifstream in;
    getTestFile("PhotoPair_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    double multiplier;
    bool lpm;
    double density_correction;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> multiplier >> lpm >> density_correction >> energy
        >> parametrization >> dNdx_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["lpm"] = lpm;

        auto cross
            = make_photopairproduction(particle_def, *medium, true, config, density_correction);

        dNdx_new = cross->CalculatedNdx(energy);

        // low-energy corrections make interpolations harder
        if (parametrization == "kochmotz" && energy < 50.)
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-2 * dNdx_stored);
        else
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
