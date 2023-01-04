#include "gtest/gtest.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/Factories/ComptonFactory.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSALTestUtilities/TestFilesHandling.h"

#include <spdlog/spdlog.h>
#include <nlohmann/json.hpp>

using namespace PROPOSAL;

TEST(Compton, Test_of_dEdx)
{
    std::ifstream in;
    getTestFile("Compton_dEdx.txt", in);
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    std::cout.precision(16);

    while (in >> mediumName >> ecut >> vcut >> multiplier >> energy
        >> dEdx_stored >> parametrization) {

        ParticleDef particle_def = GammaDef();
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_compton(particle_def, *medium, ecuts, false, config);

        dEdx_new = cross->CalculatedEdx(energy);

        ASSERT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);
    }
}

TEST(Compton, Test_of_dNdx)
{
    std::ifstream in;
    getTestFile("Compton_dNdx.txt", in);

    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    std::cout.precision(16);

    while (in >> mediumName >> ecut >> vcut >> multiplier >> energy
        >> dNdx_stored >> parametrization) {

        ParticleDef particle_def = GammaDef();
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_compton(particle_def, *medium, ecuts, false, config);

        dNdx_new = cross->CalculatedNdx(energy);

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
    }
}

TEST(Compton, Test_of_e)
{
    std::ifstream in;
    getTestFile("Compton_e.txt", in);

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

    std::cout.precision(16);

    RandomGenerator::Get().SetSeed(0);

    while (in >> mediumName >> ecut >> vcut >> multiplier >> energy >> rnd1
        >> rnd2 >> stochastic_loss_stored >> parametrization) {

        ParticleDef particle_def = GammaDef();
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_compton(particle_def, *medium, ecuts, false, config);

        auto components = medium->GetComponents();
        auto comp = components.at(int(rnd2 * components.size()));

        auto dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());

        if ( ecut == "inf" && vcut == 1 ) {
            EXPECT_THROW(cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp), std::logic_error);
        } else {
            auto stochastic_loss = cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp);
            EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-4 * stochastic_loss_stored);
            // cross check
            if (dNdx_for_comp > 0) {
                auto rate_rnd = cross->CalculateCumulativeCrosssection(energy, comp.GetHash(), stochastic_loss);
                EXPECT_NEAR(rate_rnd/dNdx_for_comp, rnd1, 1e-3);
            }
        }
    }
}

TEST(Compton, Test_of_dEdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Compton_dEdx.txt", in);

    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    std::cout.precision(16);

    while (in >> mediumName >> ecut >> vcut >> multiplier >> energy
        >> dEdx_stored >> parametrization) {

        ParticleDef particle_def = GammaDef();
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_compton(particle_def, *medium, ecuts, true, config);

        dEdx_new = cross->CalculatedEdx(energy);

        if (vcut * energy == std::stod(ecut))
            // expecting a kink here (issue #250)
            EXPECT_NEAR(dEdx_new, dEdx_stored, 5e-2 * dEdx_stored);
        else
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-4 * dEdx_stored);
    }
}

TEST(Compton, Test_of_dNdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Compton_dNdx.txt", in);

    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    std::cout.precision(16);

    while (in >> mediumName >> ecut >> vcut >> multiplier >> energy
        >> dNdx_stored >> parametrization) {

        ParticleDef particle_def = GammaDef();
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_compton(particle_def, *medium, ecuts, true, config);

        dNdx_new = cross->CalculatedNdx(energy);

        EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
    }
}

TEST(Compton, Test_of_e_Interpolant)
{
    std::ifstream in;
    getTestFile("Compton_e.txt", in);

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

    std::cout.precision(16);

    RandomGenerator::Get().SetSeed(0);

    while (in >> mediumName >> ecut >> vcut >> multiplier >> energy >> rnd1
        >> rnd2 >> stochastic_loss_stored >> parametrization) {

        ParticleDef particle_def = GammaDef();
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_compton(particle_def, *medium, ecuts, true, config);

        auto components = medium->GetComponents();
        auto comp = components.at(int(rnd2 * components.size()));

        auto dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());

        if ( ecut == "inf" && vcut == 1 ) {
            EXPECT_THROW(cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp), std::logic_error);
        } else {
            auto stochastic_loss = cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp);
            if (vcut * energy == std::stod(ecut) || rnd1 < 0.05)
                // expecting a kink here (issue #250)
                // The lower edge of the kinematic range is poorly interpolated
                // due to a discontinuity (issue #250)
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 5e-2 * stochastic_loss_stored);
            else
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-3 * stochastic_loss_stored);

            // cross check
            if (dNdx_for_comp > 0) {
                auto rate_rnd = cross->CalculateCumulativeCrosssection(energy, comp.GetHash(), stochastic_loss);
                EXPECT_NEAR(rate_rnd/dNdx_for_comp, rnd1, 1e-3);
            }
        }
    }
}

int main(int argc, char** argv)
{
    if (argc > 1) {
        auto debug = std::atof(argv[1]);
        if (debug)
            Logging::SetGlobalLoglevel(spdlog::level::level_enum::trace);
    }
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
