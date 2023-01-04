
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
    std::string ecut;
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
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

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
    std::string ecut;
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
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

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
    std::string ecut;
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
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["hard_component"] = hard_component;

        auto cross = make_photonuclearreal(
            particle_def, *medium, ecuts, false, config);

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
                if (energy >= 1e10 && hard_component == true) {
                    // Integration issues in dNdx due to a discontinuity in the
                    // differential cross section (issue #124)
                    EXPECT_NEAR(rate_rnd/dNdx_for_comp, rnd1, 5e-3);
                } else {
                    EXPECT_NEAR(rate_rnd/dNdx_for_comp, rnd1, 1e-3);
                }
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
    std::string ecut;
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
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["hard_component"] = hard_component;

        auto cross
            = make_photonuclearreal(particle_def, *medium, ecuts, true, config);

        dEdx_new = cross->CalculatedEdx(energy);

        if (hard_component == 1 && energy == 1e5) {
            // kink in differential cross section (see issue #124)
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-1 * dEdx_stored);
        } else if (vcut * energy == std::stod(ecut)) {
            // kink in interpolated function (issue #250)
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-1 * dEdx_stored);
        } else {
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);
        }
    }
}

TEST(PhotoRealPhotonAssumption, Test_of_dNdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Photo_Real_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
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
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["hard_component"] = hard_component;

        auto cross
            = make_photonuclearreal(particle_def, *medium, ecuts, true, config);

        dNdx_new = cross->CalculatedNdx(energy);
        if (energy * vcut == std::stod(ecut)) {
            // kink in interpolated function (issue #250)
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-1 * dNdx_stored);
        } else if (hard_component == 1 && energy >= 1e10) {
            // for high E, the high_component correction produces artefacts due
            // to a hard cutoff for v=1e-7 in the differential crosssection
            // (issue #124)
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-2 * dNdx_stored);
        } else if (hard_component == 1 && energy == 1e5) {
            // kink in differential crosssection when hard_component is enabled
            // at around 1e5 MeV (issue #124)
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-1 * dNdx_stored);
        } else if (parametrization == "Rhode" && energy == 1e4) {
            // Rhode parametrization hard to interpolate for this energy range
            // due to its complicated structure
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
    std::string ecut;
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
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

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

        if ( ecut == "inf" && vcut == 1 ) {
            EXPECT_THROW(cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp), std::logic_error);
        } else {
            auto stochastic_loss = cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp);

            if ( rnd1 < 0.05 && std::stod(ecut) == 500) {
                // Not enough nodes at very small losses!
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 5e-2 * stochastic_loss_stored);
            } else if (rnd1 < 0.02) {
                // The lower edge of the kinematic range is poorly interpolated
                // due to a discontinuity for v -> v_cut (issue #250)
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-2 * stochastic_loss_stored);
            } else if (hard_component == 1 && energy >= 1e10) {
                // for high E, the high_component correction produces artefacts due
                // to a hard cutoff for v=1e-7 in the differential crosssection
                // (issue #124)
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-1 * stochastic_loss_stored);
            } else if (energy >= 1e10 && std::stod(ecut) == 500) {
                // for high E, there are integration problems also without the
                // hard component, causing problems in the dNdx interpolation
                // (issue #124)
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-1 * stochastic_loss_stored);
            } else if (particleName == "TauMinus" && energy >= 1e8 && std::stod(ecut) == 500) {
                // same issue as above, but it starts earlier when taus are involved
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-2 * stochastic_loss_stored);
            } else if (hard_component == 1 && energy == 1e5) {
                    // kink in differential crosssection when hard_component is enabled
                    // at around 1e5 MeV (issue #124)
                    EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 5e-2 * stochastic_loss_stored);
            } else if (parametrization == "Rhode" && energy == 1e4) {
                // Rhode parametrization hard to interpolate for this energy range
                // due to its complicated structure
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 5e-2 * stochastic_loss_stored);
            } else if (vcut * energy == std::stod(ecut)) {
                // kink in interpolated function (issue #250)
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-2 * stochastic_loss_stored);
            } else if (rnd1 > 0.98) {
                // inaccurate dNdx interpolant for v->1
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-1 * stochastic_loss_stored);
            } else if (particleName == "EMinus" && std::stod(ecut) == 500) {
                // Jump in integrated functions are seen almost everywhere for
                // electrons (issue #124)
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 5e-3 * stochastic_loss_stored);
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
    std::string ecut;
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
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

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
    std::string ecut;
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
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

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
    std::string ecut;
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
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["shadow"] = shadowing;

        auto cross
            = make_photonuclearQ2(particle_def, *medium, ecuts, false, config);

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

TEST(PhotoQ2Integration, Test_of_dEdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Photo_Q2_dEdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
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
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["shadow"] = shadowing;

        auto cross
            = make_photonuclearQ2(particle_def, *medium, ecuts, true, config);

        dEdx_new = cross->CalculatedEdx(energy);
        if (std::stod(ecut) == vcut * energy) {
            // kink in interpolated function (issue #250)
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-1 * dEdx_stored);
        } else {
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);
        }
    }
}

TEST(PhotoQ2Integration, Test_of_dNdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Photo_Q2_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
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
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["shadow"] = shadowing;

        auto cross
            = make_photonuclearQ2(particle_def, *medium, ecuts, true, config);

        dNdx_new = cross->CalculatedNdx(energy);
        if (energy * vcut == std::stod(ecut)) {
            // kink in interpolated function (issue #250)
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-1 * dNdx_stored);
        } else {
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
        }
    }
}

TEST(PhotoQ2Integration, Test_of_e_Interpolant)
{
    std::ifstream in;
    getTestFile("Photo_Q2_e.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
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
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["shadow"] = shadowing;

        auto cross
            = make_photonuclearQ2(particle_def, *medium, ecuts, true, config);

        auto components = medium->GetComponents();

        auto comp = components.at(int(rnd2 * components.size()));

        auto dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());

        if ( ecut == "inf" && vcut == 1 ) {
            EXPECT_THROW(cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp), std::logic_error);
        } else {
            auto stochastic_loss = cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp);

            if (energy * vcut == std::stod(ecut)) {
                // kink in interpolated function (issue #250)
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-1 * stochastic_loss_stored);
            } else if (rnd1 < 0.07 || rnd1 > 0.995) {
                // The lower edge of the kinematic range is poorly interpolated
                // due to a discontinuity (issue #250)
                // inaccurate dNdx interpolant for v->1
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 5e-2 * stochastic_loss_stored);
            } else if (parametrization == "AbtFT" && energy >= 1e10) {
                // problems with integration, leading to kinks in integrated
                // function for high energies (issue #124)
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-2 * stochastic_loss_stored);
            } else if (parametrization == "ButkevichMikheyev" && energy == 100000000000 && rnd1 < 0.1) {
                // combination of issue #250 and issue #124
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 5e-3 * stochastic_loss_stored);
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
