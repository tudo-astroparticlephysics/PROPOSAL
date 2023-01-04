
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/Factories/IonizationFactory.h"
#include "PROPOSALTestUtilities/TestFilesHandling.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "MuMinus")
    {
        return MuMinusDef();
    } else if (name == "TauMinus")
    {
        return TauMinusDef();
    } else
    {
        return EMinusDef();
    }
}

TEST(Ionization, Test_of_dEdx)
{
    std::ifstream in;
    getTestFile("Ioniz_dEdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> dEdx_stored >> parametrization)
    {
        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_ionization(particle_def, *medium, ecuts, false,
                                     config);

        dEdx_new = cross->CalculatedEdx(energy);
        EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-6 * dEdx_stored);
    }
}

TEST(Ionization, Test_of_dNdx)
{
    std::ifstream in;
    getTestFile("Ioniz_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> dNdx_stored >> parametrization)
    {
        ParticleDef particle_def = getParticleDef(particleName);

        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_ionization(particle_def, *medium, ecuts, false,
                                     config);

        dNdx_new = cross->CalculatedNdx(energy);
        EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-6 * dNdx_stored);
    }
}

TEST(Ionization, Test_Stochastic_Loss)
{
    std::ifstream in;
    getTestFile("Ioniz_e.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;

    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> rnd1 >> rnd2 >> stochastic_loss_stored >> parametrization) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_ionization(particle_def, *medium, ecuts, false,
                                     config);

        auto dNdx = cross->CalculatedNdx(energy);

        if ( ecut == "inf" && vcut == 1 ) {
            EXPECT_THROW(cross->CalculateStochasticLoss(medium->GetHash(), energy, rnd1 * dNdx), std::logic_error);
        } else {
            auto stochastic_loss = cross->CalculateStochasticLoss(medium->GetHash(), energy, rnd1 * dNdx);
            EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-6 * stochastic_loss_stored);

            // cross check
            if (dNdx > 0) {
                auto rate_rnd = cross->CalculateCumulativeCrosssection(energy, medium->GetHash(), stochastic_loss);
                EXPECT_NEAR(rate_rnd/dNdx, rnd1, 1e-3);
            }
        }
    }
}

TEST(Ionization, Test_of_dEdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Ioniz_dEdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> dEdx_stored >> parametrization) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_ionization(particle_def, *medium, ecuts, true,
                                     config);

        dEdx_new = cross->CalculatedEdx(energy);
        if (vcut * energy == std::stod(ecut)) {
            // kink in interpolated function (issue #250)
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-2 * dEdx_stored);
        } else {
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);
        }

    }
}

TEST(Ionization, Test_of_dNdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Ioniz_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> dNdx_stored >> parametrization) {

        ParticleDef particle_def = getParticleDef(particleName);

        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_ionization(particle_def, *medium, ecuts, true,
                                     config);

        dNdx_new = cross->CalculatedNdx(energy);
        if (vcut * energy == std::stod(ecut)) {
            // kink in interpolated function (issue #250)
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-1 * dNdx_stored);
        } else {
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
        }
    }
}

TEST(Ionization, Test_of_e_Interpolant)
{
    std::ifstream in;
    getTestFile("Ioniz_e.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;

    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> rnd1 >> rnd2 >> stochastic_loss_stored >> parametrization) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_ionization(particle_def, *medium, ecuts, true,
                                     config);

        auto dNdx = cross->CalculatedNdx(energy);

        if ( ecut == "inf" && vcut == 1 || dNdx == 0) {
            EXPECT_THROW(cross->CalculateStochasticLoss(medium->GetHash(), energy, rnd1 * dNdx), std::logic_error);
        } else {
            auto stochastic_loss = cross->CalculateStochasticLoss(medium->GetHash(), energy, rnd1 * dNdx);
            if (vcut * energy == std::stod(ecut)) {
                // kink in interpolated function (issue #250)
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 5e-2 * stochastic_loss_stored);
            } else if (rnd1 > 0.99) {
                // Integration routine unstable here for very high rnd numbers (issue #254)
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-1 * stochastic_loss_stored);
            } else {
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-3 * stochastic_loss_stored);
            }

            // cross check
            if (dNdx > 0) {
                auto rate_rnd = cross->CalculateCumulativeCrosssection(energy, medium->GetHash(), stochastic_loss);
                EXPECT_NEAR(rate_rnd/dNdx, rnd1, 1e-3);
            }
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
