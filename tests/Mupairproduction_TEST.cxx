
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/Factories/MupairProductionFactory.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/secondaries/parametrization/mupairproduction/KelnerKokoulinPetrukhinMupairProduction.h"
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

TEST(Mupairproduction, Test_of_dEdx)
{
    std::ifstream in;
    getTestFile("Mupair_dEdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    double energy;
    std::string parametrization;
    double dEdx_stored;
    double dEdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> parametrization >> dEdx_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_mupairproduction(particle_def, *medium, ecuts, false,
                                           config);

        dEdx_new = cross->CalculatedEdx(energy);
        EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-6 * dEdx_stored);
    }
}

TEST(Mupairproduction, Test_of_dNdx)
{
    std::ifstream in;
    getTestFile("Mupair_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> parametrization >> dNdx_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_mupairproduction(particle_def, *medium, ecuts, false,
                                           config);

        dNdx_new = cross->CalculatedNdx(energy);
        EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-6 * dNdx_stored);
    }
}

TEST(Mupairproduction, Test_Stochastic_Loss)
{
    std::ifstream in;
    getTestFile("Mupair_e.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    double energy;
    std::string parametrization;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> parametrization >> rnd1 >> rnd2 >> stochastic_loss_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_mupairproduction(particle_def, *medium, ecuts, false,
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

TEST(Mupairproduction, Test_Calculate_Rho)
{
    std::ifstream in;
    getTestFile("Mupair_rho.txt", in);

    std::string particleName;
    std::string mediumName;
    double v;
    double energy;
    std::string parametrization;
    double rnd1, rnd2;
    double E1_stored;
    double E2_stored;
    double E1_new;
    double E2_new;
    double rho;

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> v >> energy >> rnd1 >> rnd2
        >> E1_stored >> E2_stored) {

    ParticleDef particle_def = getParticleDef(particleName);
    std::shared_ptr<const Medium> medium = CreateMedium(mediumName);

    auto fac = secondaries::KelnerKokoulinPetrukhinMupairProduction(particle_def, *medium);
    rho = fac.CalculateRho(energy, v, medium->GetComponents().front(), rnd1, rnd2);

    E1_new = 0.5 * v * energy * (1 + rho);
    E2_new = 0.5 * v * energy * (1 - rho);

    EXPECT_NEAR(E1_new, E1_stored, 1E-3 * E1_stored);
    EXPECT_NEAR(E2_new, E2_stored, 1E-3 * E2_stored);
    }
}

TEST(Mupairproduction, Test_of_dEdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Mupair_dEdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    double energy;
    std::string parametrization;
    double dEdx_stored;
    double dEdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut
        >> multiplier >> energy >> parametrization >> dEdx_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_mupairproduction(particle_def, *medium, ecuts, true,
                                           config);

        dEdx_new = cross->CalculatedEdx(energy);

        if (energy * vcut == std::stod(ecut)) {
            // kink in interpolated function (issue #250)
            EXPECT_NEAR(dEdx_new, dEdx_stored, 2e-1 * dEdx_stored);
        } else {
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);
        }
    }
}

TEST(Mupairproduction, Test_of_dNdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Mupair_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> parametrization >> dNdx_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_mupairproduction(particle_def, *medium, ecuts, true,
                                           config);
        dNdx_new = cross->CalculatedNdx(energy);

        if (energy * vcut == std::stod(ecut)) {
            // kink in interpolated function (issue #250)
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-1 * dNdx_stored);
        } else {
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
        }
    }
}

TEST(Mupairproduction, Test_of_e_Interpolant)
{
    std::ifstream in;
    getTestFile("Mupair_e.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    double energy;
    std::string parametrization;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;

    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> parametrization >> rnd1 >> rnd2 >> stochastic_loss_stored) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_mupairproduction(particle_def, *medium, ecuts, true,
                                           config);

        auto components = medium->GetComponents();
        auto comp = components.at(int(rnd2 * components.size()));

        auto dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());

        if ( ecut == "inf" && vcut == 1 ) {
            EXPECT_THROW(cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp), std::logic_error);
        } else {
            auto stochastic_loss = cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp);
            if (energy * vcut == std::stod(ecut)) {
                // kink in interpolated function (issue #250)
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 1e-2 * stochastic_loss_stored);
            } else if (rnd1 > 0.95) {
                // dNdx interpolation for v->1 is not accurate enough (issue #253)
                EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, 5e-2 * stochastic_loss_stored);
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
