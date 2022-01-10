
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSALTestUtilities/TestFilesHandling.h"
#include "PROPOSAL/crosssection/Factories/PhotonuclearFactory.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include <fstream>

using namespace PROPOSAL;

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "MuMinus") {
        return MuMinusDef();
    } else if (name == "TauMinus") {
        return TauMinusDef();
    } else {
        return EMinusDef();
    }
}

TEST(PhotoRealPhotonAssumption, Test_of_dEdx)
{
    std::ifstream in;
    getTestFile("Photo_Real_dEdx.txt", in);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    bool hard_component;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> dEdx_stored >> parametrization >> hard_component) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["hard_component"] = hard_component;

        auto cross = make_photonuclearreal(
            particle_def, *medium, ecuts, false, config);

        dEdx_new = cross->CalculatedEdx(energy);

        EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-5 * dEdx_stored);
    }
}

TEST(PhotoRealPhotonAssumption, Test_of_dNdx)
{
    std::ifstream in;
    getTestFile("Photo_Real_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    bool hard_component;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> dNdx_stored >> parametrization >> hard_component) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["hard_component"] = hard_component;

        auto cross = make_photonuclearreal(
            particle_def, *medium, ecuts, false, config);

        dNdx_new = cross->CalculatedNdx(energy);

        EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-5 * dNdx_stored);
    }
}

TEST(PhotoRealPhotonAssumption, Test_of_e)
{
    std::ifstream in;
    getTestFile("Photo_Real_e.txt", in);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    bool hard_component;
    double energy;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> rnd1 >> rnd2 >> stochastic_loss_stored >> parametrization
        >> hard_component) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["hard_component"] = hard_component;

        auto cross = make_photonuclearreal(
            particle_def, *medium, ecuts, false, config);

        auto components = medium->GetComponents();

        auto comp = components.at(int(rnd2 * components.size()));

        auto dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());

        if ( ecut == INF && vcut == 1 ) {
            EXPECT_THROW(cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp), std::logic_error);
        } else {
            auto stochastic_loss = cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp);
            EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-3 * stochastic_loss_stored);

            // cross check
            if (dNdx_for_comp > 0) {
                auto rate_rnd = cross->CalculateCumulativeCrosssection(energy, comp.GetHash(), stochastic_loss);
                EXPECT_NEAR(rate_rnd/dNdx_for_comp, rnd1, 5e-3);
            }
        }
    }
}

TEST(PhotoRealPhotonAssumption, Test_of_dEdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Photo_Real_dEdx.txt", in);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    bool hard_component;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> dEdx_stored >> parametrization >> hard_component) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["hard_component"] = hard_component;

        auto cross
            = make_photonuclearreal(particle_def, *medium, ecuts, true, config);

        dEdx_new = cross->CalculatedEdx(energy);

        if (hard_component == 1 && energy == 1e5)
            EXPECT_NEAR(dEdx_new, dEdx_stored,
                1e-1 * dEdx_stored); // kink in function (see issue #124)
        else if (vcut * energy == ecut)
            EXPECT_NEAR(
                dEdx_new, dEdx_stored, 1e-1 * dEdx_stored); // kink in function
        else if (parametrization == "Rhode" && energy <= 10000)
            EXPECT_NEAR(dEdx_new, dEdx_stored,
                1e-2 * dEdx_stored); // slight "bumps" in function, hard to
                                     // interpolate
        else
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);
    }
}

TEST(PhotoRealPhotonAssumption, Test_of_dNdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Photo_Real_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    bool hard_component;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> dNdx_stored >> parametrization >> hard_component) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["hard_component"] = hard_component;

        auto cross
            = make_photonuclearreal(particle_def, *medium, ecuts, true, config);

        dNdx_new = cross->CalculatedNdx(energy);
        if (energy * vcut == ecut) {
            EXPECT_NEAR(
                    dNdx_new, dNdx_stored, 1e-1 * dNdx_stored); // kink in function
        } else if (hard_component == 1 && energy >= 1e10) {
            EXPECT_NEAR(dNdx_new, dNdx_stored,
                        1e-2 * dNdx_stored); // for high E, high_component correction
                        // produces artefacts due to hard cutoff
                        // for v=1e-7 in differential cross section
                        // (see issue #124)
        } else if (hard_component == 1 && energy == 1e5) {
            EXPECT_NEAR(dNdx_new, dNdx_stored,
                        1e-1 * dNdx_stored); // kink in hard_component for E=1e5 (see
                        // issue #124)
        } else if (parametrization == "Rhode" && energy == 1e4) {
            // rhode parametrization hard to interpolate for this energy range
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-2 * dNdx_stored);
        } else {
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
        }
    }
}

TEST(PhotoRealPhotonAssumption, Test_of_e_Interpolant)
{
    std::ifstream in;
    getTestFile("Photo_Real_e.txt", in);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    bool hard_component;
    double energy;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> rnd1 >> rnd2 >> stochastic_loss_stored >> parametrization
        >> hard_component) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["hard_component"] = hard_component;

        auto cross
            = make_photonuclearreal(particle_def, *medium, ecuts, true, config);
        auto cross_integral
            = make_photonuclearreal(particle_def, *medium, ecuts, true, config);

        auto components = medium->GetComponents();

        auto comp = components.at(int(rnd2 * components.size()));

        auto dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());

        if ( ecut == INF && vcut == 1 ) {
            EXPECT_THROW(cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp), std::logic_error);
        } else {
            auto stochastic_loss = cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp);
            if ( rnd1 < 0.05 || rnd1 > 0.95) {
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-2 * stochastic_loss_stored);
            } else if (hard_component == 1 && energy >= 1e10) {
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-2 * stochastic_loss_stored);
            } else if (hard_component == 1 && energy == 1e5) {
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-2 * stochastic_loss_stored);
            } else if (parametrization == "Rhode" && energy == 1e4) {
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-2 * stochastic_loss_stored);
            } else if (vcut * energy == ecut) {
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

// Photonuclear parametrization with Q2 Integration

TEST(PhotoQ2Integration, Test_of_dEdx)
{
    std::ifstream in;
    getTestFile("Photo_Q2_dEdx.txt", in);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    std::string shadowing;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> dEdx_stored >> parametrization >> shadowing) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["shadow"] = shadowing;

        auto cross
            = make_photonuclearQ2(particle_def, *medium, ecuts, false, config);

        dEdx_new = cross->CalculatedEdx(energy);

        EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-5 * dEdx_stored);
    }
}

TEST(PhotoQ2Integration, Test_of_dNdx)
{
    std::ifstream in;
    getTestFile("Photo_Q2_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    std::string shadowing;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> dNdx_stored >> parametrization >> shadowing) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["shadow"] = shadowing;

        auto cross
            = make_photonuclearQ2(particle_def, *medium, ecuts, false, config);

        dNdx_new = cross->CalculatedNdx(energy);

        EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-5 * dNdx_stored);
    }
}

TEST(PhotoQ2Integration, Test_of_e)
{
    std::ifstream in;
    getTestFile("Photo_Q2_e.txt", in);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    std::string shadowing;
    double energy;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> rnd1 >> rnd2 >> stochastic_loss_stored >> parametrization
        >> shadowing) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["shadow"] = shadowing;

        auto cross
            = make_photonuclearQ2(particle_def, *medium, ecuts, false, config);

        auto components = medium->GetComponents();

        auto comp = components.at(int(rnd2 * components.size()));

        auto dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());

        if ( ecut == INF && vcut == 1 ) {
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

TEST(PhotoQ2Integration, Test_of_dEdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Photo_Q2_dEdx.txt", in);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    std::string shadowing;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> dEdx_stored >> parametrization >> shadowing) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["shadow"] = shadowing;

        auto cross
            = make_photonuclearQ2(particle_def, *medium, ecuts, true, config);

        dEdx_new = cross->CalculatedEdx(energy);
        if (ecut == vcut * energy)
            EXPECT_NEAR(
                dEdx_new, dEdx_stored, 1e-1 * dEdx_stored); // kink in integral
        else if (parametrization == "AbramowiczLevinLevyMaor91"
            && shadowing == "ButkevichMikheyev" && particleName == "MuMinus"
            && energy == 100000000000 && ecut == INF && vcut == 1)
            EXPECT_NEAR(dEdx_new, dEdx_stored,
                1e-2 * dEdx_stored); // bump in integral for these specific
                                     // parameter
        else
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);
    }
}

TEST(PhotoQ2Integration, Test_of_dNdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Photo_Q2_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    std::string shadowing;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> dNdx_stored >> parametrization >> shadowing) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["shadow"] = shadowing;

        auto cross
            = make_photonuclearQ2(particle_def, *medium, ecuts, true, config);

        dNdx_new = cross->CalculatedNdx(energy);
        if (energy * vcut == ecut)
            EXPECT_NEAR(
                dNdx_new, dNdx_stored, 1e-1 * dNdx_stored); // kink in integral
        else
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
    }
}

TEST(PhotoQ2Integration, Test_of_e_Interpolant)
{
    std::ifstream in;
    getTestFile("Photo_Q2_e.txt", in);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    std::string shadowing;
    double energy;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> rnd1 >> rnd2 >> stochastic_loss_stored >> parametrization
        >> shadowing) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["shadow"] = shadowing;

        auto cross
            = make_photonuclearQ2(particle_def, *medium, ecuts, true, config);

        auto components = medium->GetComponents();

        auto comp = components.at(int(rnd2 * components.size()));

        auto dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());

        if ( ecut == INF && vcut == 1 ) {
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

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
