
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/Factories/EpairProductionFactory.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSALTestUtilities/TestFilesHandling.h"

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

TEST(Epairproduction, Test_of_dEdx)
{
    std::ifstream in;
    getTestFile("Epair_dEdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double dEdx_stored;
    double dEdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> lpm >> energy >> parametrization >> dEdx_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["lpm"] = lpm;

        auto cross = make_epairproduction(particle_def, *medium, ecuts,
                                          false, config);

        dEdx_new = cross->CalculatedEdx(energy);
        EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-5 * dEdx_stored);
    }
}

TEST(Epairproduction, Test_of_dNdx)
{
    std::ifstream in;
    getTestFile("Epair_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >>
        lpm >> energy >> parametrization >> dNdx_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["lpm"] = lpm;

        auto cross = make_epairproduction(particle_def, *medium, ecuts, false,
                                          config);

        dNdx_new = cross->CalculatedNdx(energy);

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);
    }
}

TEST(Epairproduction, Test_Stochastic_Loss)
{
    std::ifstream in;
    getTestFile("Epair_e.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >>
        lpm >> energy >> parametrization >> rnd1 >> rnd2 >> stochastic_loss_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["lpm"] = lpm;

        auto cross = make_epairproduction(particle_def, *medium, ecuts, false,
                                          config);

        auto components = medium->GetComponents();
        auto comp = components.at(int(rnd2 * components.size()));

        auto dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());

        if ( ecut == "inf" && vcut == 1 ) {
            EXPECT_THROW(cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp), std::logic_error);
        } else {
            auto stochastic_loss = cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp);
            EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-3 * stochastic_loss_stored);

            // cross check
            if (dNdx_for_comp > 0) {
                auto rate_rnd = cross->CalculateCumulativeCrosssection(energy, comp.GetHash(), stochastic_loss);
                EXPECT_NEAR(rate_rnd/dNdx_for_comp, rnd1, 1e-3);
            }
        }
    }
}

TEST(Epairproduction, Test_of_dEdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Epair_dEdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double dEdx_stored;
    double dEdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> lpm >> energy >> parametrization >> dEdx_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["lpm"] = lpm;

        auto cross = make_epairproduction(particle_def, *medium, ecuts, true,
                                          config);

        dEdx_new = cross->CalculatedEdx(energy);

        if (particleName == "TauMinus" && mediumName == "uranium" && energy == 1e4) {
            // Transition from zero to non-zero values is not continuous, causing
            // interpolation problems in uranium around 1e4 MeV (issue #250)
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e0 * dEdx_stored);
        } else if (vcut * energy == std::stod(ecut)) {
            // kink in interpolated function (issue #250)
            EXPECT_NEAR(dEdx_new, dEdx_stored, 5e-2 * dEdx_stored);
        } else {
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);
        }
    }
}

TEST(Epairproduction, Test_of_dNdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Epair_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> lpm >> energy >> parametrization >> dNdx_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["lpm"] = lpm;

        auto cross = make_epairproduction(particle_def, *medium, ecuts, true,
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

TEST(Epairproduction, Test_of_e_Interpolant)
{
    std::ifstream in;
    getTestFile("Epair_e.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;

    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> lpm >> energy >> parametrization >> rnd1 >> rnd2 >> stochastic_loss_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["lpm"] = lpm;

        auto cross = make_epairproduction(particle_def, *medium, ecuts, true,
                                          config);

        auto components = medium->GetComponents();
        auto comp = components.at(int(rnd2 * components.size()));

        auto dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());

        if ( ecut == "inf" && vcut == 1 || dNdx_for_comp == 0 ) {
            EXPECT_THROW(cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp), std::logic_error);
        } else {
            auto stochastic_loss = cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp);

            if (energy * vcut == std::stod(ecut)) {
                // kink in interpolated function (issue #250)
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 5e-2 * stochastic_loss_stored);
            } else if (rnd1 < 0.1) {
                // The lower edge of the kinematic range is poorly interpolated
                // due to a discontinuity (issue #250)
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 5e-2 * stochastic_loss_stored);
            } else if (particleName == "TauMinus" && rnd1 < 0.2) {
                // Same as above, but the range is a bit higher for tau leptons
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 5e-2 * stochastic_loss_stored);
            } else if (particleName == "TauMinus" && energy <= 1e5 && rnd1 > 0.9){
                // The transition from zero to non-zero values in dNdx is not continuous,
                // causing the dNdx interpolation to be worse than usual (issue #250)
                // This error is amplified for high stochastic losses because the losses
                // are very sensitive to changes in the cumulative crosssection.
                // This causes the losses to be less accurate than required.
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-2 * stochastic_loss_stored);
            } else if (particleName == "EMinus" && energy >= 1e10) {
                // Artefacts in dNdx integral, causing the dNdx interpolation to
                // have an uncertainty, which causes the stochastic losses to be
                // not accurate enough.
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-2 * stochastic_loss_stored);
            } else {
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-3 * stochastic_loss_stored);
            }

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
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
