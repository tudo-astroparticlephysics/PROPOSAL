#include "gtest/gtest.h"
#include <nlohmann/json.hpp>

#include "PROPOSALTestUtilities/TestFilesHandling.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/Factories/BremsstrahlungFactory.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

constexpr static double interpolation_precision = 1.e-3;

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

TEST(Bremsstrahlung, Test_of_dEdx)
{
    std::ifstream in;
    getTestFile("Brems_dEdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    std::string parametrization;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> dEdx_stored >>
        parametrization) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["lpm"] = lpm;

        auto cross = make_bremsstrahlung(particle_def, *medium, ecuts, false,
                                         config);

        dEdx_new = cross->CalculatedEdx(energy);

        EXPECT_NEAR(dEdx_new, dEdx_stored, interpolation_precision * dEdx_stored);
    }
}

TEST(Bremsstrahlung, Test_of_dNdx)
{
    std::ifstream in;
    getTestFile("Brems_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    std::string parametrization;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> dNdx_stored >>
        parametrization) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["lpm"] = lpm;

        auto cross = make_bremsstrahlung(particle_def, *medium, ecuts, false,
                                         config);

        dNdx_new = cross->CalculatedNdx(energy);

        EXPECT_NEAR(dNdx_new, dNdx_stored, interpolation_precision * dNdx_stored);
    }
}

TEST(Bremsstrahlung, Test_of_e)
{
    std::ifstream in;
    getTestFile("Brems_e.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    std::string parametrization;
    double energy;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;

    std::cout.precision(16);

    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> rnd1 >> rnd2 >>
        stochastic_loss_stored >> parametrization) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["lpm"] = lpm;

        auto cross = make_bremsstrahlung(particle_def, *medium, ecuts, false,
                                         config);

        auto components = medium->GetComponents();
        auto comp = components.at(int(rnd2 * components.size()));

        auto dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());

        if ( ecut == "inf" && vcut == 1 ) {
            EXPECT_THROW(cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp), std::logic_error);
        } else {
            auto stochastic_loss = cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp);
            EXPECT_NEAR(stochastic_loss, stochastic_loss_stored, interpolation_precision*stochastic_loss_stored);

            // cross check
            if (dNdx_for_comp > 0) {
                auto rate_rnd = cross->CalculateCumulativeCrosssection(energy, comp.GetHash(), stochastic_loss);
                // for high-energy electrons with lpm enabled, the differential crosssetion rises for v->1, so the
                // differential cross section becomes harder to integrate (issue #123)
                if (particleName == "EMinus" && energy >= 1e10 && lpm == true)
                    EXPECT_NEAR(rate_rnd/dNdx_for_comp, rnd1, 5e-2);
                else
                    EXPECT_NEAR(rate_rnd/dNdx_for_comp, rnd1, 1e-3);
            }
        }


    }
}

TEST(Bremsstrahlung, Test_of_dEdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Brems_dEdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    std::string parametrization;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >>
    dEdx_stored >> parametrization) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["lpm"] = lpm;

        auto cross = make_bremsstrahlung(particle_def, *medium, ecuts, true,
                                         config);

        dEdx_new = cross->CalculatedEdx(energy);


        if (particleName == "TauMinus" && energy < 1e5) {
            // For taus, the kinematic upper limit (v_max) introduces kinks in
            // the function that we need to integrate. However, bremsstrahlung
            // effects for taus are negligible for these energies (issue #250)
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e0 * dEdx_stored);
        } else if (vcut * energy == std::stod(ecut)) {
            // expecting a kink here (issue #250)
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-1 * dEdx_stored);
        } else if (particleName == "EMinus" && mediumName == "uranium" && energy == 1e10 && lpm == true) {
            // There is a small discontinuity in the function that
            // is hard to interpolate at exactly this energy (issue #250)
            EXPECT_NEAR(dEdx_new, dEdx_stored, 5e-3 * dEdx_stored);
        } else {
            EXPECT_NEAR(dEdx_new, dEdx_stored, interpolation_precision * dEdx_stored);
        }
    }
}

TEST(Bremsstrahlung, Test_of_dNdx_Interpolant)
{
    std::ifstream in;
    getTestFile("Brems_dNdx.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    std::string parametrization;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    std::cout.precision(16);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >>
        lpm >> energy >> dNdx_stored >> parametrization) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["lpm"] = lpm;

        auto cross = make_bremsstrahlung(particle_def, *medium, ecuts, true,
                                         config);

        dNdx_new = cross->CalculatedNdx(energy);

        if (vcut * energy == std::stod(ecut)) {
            // expecting a kink here (issue #250)
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-1 * dNdx_stored);
        } else if (particleName == "TauMinus" && energy < 1e5) {
            // For taus, the kinematic upper limit (v_max) introduces kinks in
            // the function that we need to integrate. However, bremsstrahlung
            // effects for taus are negligible for these energies (issue #250)
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-2 * dNdx_stored);
        } else {
            EXPECT_NEAR(dNdx_new, dNdx_stored, interpolation_precision * dNdx_stored);
        }
    }
}

TEST(Bremsstrahlung, Test_of_e_Interpolant)
{
    std::ifstream in;
    getTestFile("Brems_e.txt", in);

    std::string particleName;
    std::string mediumName;
    std::string ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    std::string parametrization;
    double energy;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;

    std::cout.precision(16);

    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >>
        energy >> rnd1 >> rnd2 >> stochastic_loss_stored >> parametrization) {

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(std::stod(ecut), vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["lpm"] = lpm;

        auto cross = make_bremsstrahlung(particle_def, *medium, ecuts, true,
                                         config);

        auto components = medium->GetComponents();
        auto comp = components.at(int(rnd2 * components.size()));
        auto dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());

        if (particleName == "TauMinus" && mediumName == "uranium" && energy == 1e4)
            continue; // dNdx is non-zero for the integral, but zero for the interpolant here

        if ( ecut == "inf" && vcut == 1 ) {
            EXPECT_THROW(cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp), std::logic_error);
        } else {
            auto v = cross->CalculateStochasticLoss(comp.GetHash(), energy, rnd1 * dNdx_for_comp);
            if (energy * vcut == std::stod(ecut)) {
                // expecting a kink here (issue #250)
                EXPECT_NEAR(v, stochastic_loss_stored, 1e-1 * stochastic_loss_stored);
            } else if (particleName == "TauMinus" && energy < 1e5) {
                // For taus, the kinematic upper limit (v_max) introduces kinks in
                // the function that we need to integrate. However, bremsstrahlung
                // effects for taus are negligible for these energies (issue #250)
                EXPECT_NEAR(v, stochastic_loss_stored, 1e-2 * stochastic_loss_stored);
            } else if (particleName == "EMinus" && energy >= 1e10 && lpm == true) {
                // for high-energy electrons with LPM enabled, the differential crosssetion rises for v->1, so the
                // differential cross section becomes harder to integrate (issue #123)
                EXPECT_NEAR(v, stochastic_loss_stored, 2e-1 * stochastic_loss_stored);
            } else if (particleName == "EMinus" && mediumName == "uranium" && lpm == true) {
                // Same as above, although due to the high Z of uranium, the issues
                // becomes relevant already at lower energies (issue #123)
                EXPECT_NEAR(v, stochastic_loss_stored, 1e-2 * stochastic_loss_stored);
            } else if (rnd1 > 0.93 || rnd1 < 0.04) {
                // The very high edge of the kinematic range is only poorly
                // interpolated (see issue #253)
                // The lower edge of the kinematic range is poorly interpolated
                // due to a discontinuity (issue #250)
                EXPECT_NEAR(v, stochastic_loss_stored, 2e-2 * stochastic_loss_stored);
            } else {
                EXPECT_NEAR(v, stochastic_loss_stored, interpolation_precision * stochastic_loss_stored);
            }

            // cross check (this is actually the only test we are really interested in)
            if (dNdx_for_comp > 0) {
                auto rate_rnd = cross->CalculateCumulativeCrosssection(energy, comp.GetHash(), v);
                EXPECT_NEAR(rate_rnd/dNdx_for_comp, rnd1, 1e-4);
            }

        }

    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
