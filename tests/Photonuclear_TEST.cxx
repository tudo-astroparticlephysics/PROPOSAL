
#include "gtest/gtest.h"

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

const std::string testfile_dir = "tests/TestFiles/";

// TEST(Comparison, Comparison_equal)
// {
//     ParticleDef particle_def = MuMinusDef::Get();
//     auto medium = std::make_shared<const Water>();
//     EnergyCutSettings ecuts;
//     double multiplier   = 1.;
//     bool hard_component = true;
//     ShadowButkevichMikheyev shadow;
//     InterpolationDef InterpolDef;

//     PhotoKokoulin* PhotoReal_A   = new PhotoKokoulin(particle_def, medium,
//     ecuts, multiplier, hard_component); Parametrization* PhotoReal_B = new
//     PhotoKokoulin(particle_def, medium, ecuts, multiplier, hard_component);
//     EXPECT_TRUE(*PhotoReal_A == *PhotoReal_B);

//     PhotoKokoulin param_PhotoReal(particle_def, medium, ecuts, multiplier,
//     hard_component); EXPECT_TRUE(param_PhotoReal == *PhotoReal_A);

//     PhotoIntegral* Int_PhotoReal_A        = new
//     PhotoIntegral(param_PhotoReal); CrossSectionIntegral* Int_PhotoReal_B =
//     new PhotoIntegral(param_PhotoReal); EXPECT_TRUE(*Int_PhotoReal_A ==
//     *Int_PhotoReal_B);

//     PhotoInterpolant* Interpol_PhotoReal_A        = new
//     PhotoInterpolant(param_PhotoReal, InterpolDef); CrossSectionInterpolant*
//     Interpol_PhotoReal_B = new PhotoInterpolant(param_PhotoReal,
//     InterpolDef); EXPECT_TRUE(*Interpol_PhotoReal_A ==
//     *Interpol_PhotoReal_B);

//     delete PhotoReal_A;
//     delete PhotoReal_B;
//     delete Int_PhotoReal_A;
//     delete Int_PhotoReal_B;
//     delete Interpol_PhotoReal_A;
//     delete Interpol_PhotoReal_B;

//     PhotoAbramowiczLevinLevyMaor97* PhotoQ2_A =
//         new PhotoAbramowiczLevinLevyMaor97(particle_def, medium, ecuts,
//         multiplier, shadow);
//     Parametrization* PhotoQ2_B = new
//     PhotoAbramowiczLevinLevyMaor97(particle_def, medium, ecuts, multiplier,
//     shadow); EXPECT_TRUE(*PhotoQ2_A == *PhotoQ2_B);

//     PhotoAbramowiczLevinLevyMaor97 param_Q2_integral(particle_def, medium,
//     ecuts, multiplier, shadow);
//     PhotoQ2Interpolant<PhotoAbramowiczLevinLevyMaor97>
//     param_Q2_interpol(particle_def, medium, ecuts, multiplier, shadow,
//     InterpolDef); EXPECT_TRUE(param_Q2_integral == *PhotoQ2_A);

//     PhotoIntegral* Int_PhotoQ2_A        = new
//     PhotoIntegral(param_Q2_integral); CrossSectionIntegral* Int_PhotoQ2_B =
//     new PhotoIntegral(param_Q2_integral); EXPECT_TRUE(*Int_PhotoQ2_A ==
//     *Int_PhotoQ2_B);

//     PhotoInterpolant* Interpol_PhotoQ2_A        = new
//     PhotoInterpolant(param_Q2_interpol, InterpolDef);
//     CrossSectionInterpolant* Interpol_PhotoQ2_B = new
//     PhotoInterpolant(param_Q2_interpol, InterpolDef);
//     EXPECT_TRUE(*Interpol_PhotoQ2_A == *Interpol_PhotoQ2_B);

//     delete PhotoQ2_A;
//     delete PhotoQ2_B;
//     delete Int_PhotoQ2_A;
//     delete Int_PhotoQ2_B;
//     delete Interpol_PhotoQ2_A;
//     delete Interpol_PhotoQ2_B;
// }

// TEST(Comparison, Comparison_not_equal)
// {
//     ParticleDef mu_def  = MuMinusDef::Get();
//     ParticleDef tau_def = TauMinusDef::Get();
//     auto medium_1 = std::make_shared<const Water>();
//     auto medium_2 = std::make_shared<const Ice>();
//     EnergyCutSettings ecuts_1(500, 0.05);
//     EnergyCutSettings ecuts_2(-1, 0.05);
//     double multiplier_1 = 1.;
//     double multiplier_2 = 2.;
//     bool hard_component = true;
//     ShadowButkevichMikheyev shadow_1;
//     ShadowDuttaRenoSarcevicSeckel shadow_2;
//     InterpolationDef InterpolDef;

//     PhotoKokoulin PhotoReal_A(mu_def, medium_1, ecuts_1, multiplier_1,
//     hard_component); PhotoKokoulin PhotoReal_B(tau_def, medium_1, ecuts_1,
//     multiplier_1, hard_component); PhotoKokoulin PhotoReal_C(mu_def,
//     medium_2, ecuts_1, multiplier_1, hard_component); PhotoKokoulin
//     PhotoReal_D(mu_def, medium_1, ecuts_2, multiplier_1, hard_component);
//     PhotoKokoulin PhotoReal_E(mu_def, medium_1, ecuts_1, multiplier_1,
//     !hard_component); PhotoKokoulin PhotoReal_F(mu_def, medium_1, ecuts_1,
//     multiplier_2, hard_component); EXPECT_TRUE(PhotoReal_A != PhotoReal_B);
//     EXPECT_TRUE(PhotoReal_A != PhotoReal_C);
//     EXPECT_TRUE(PhotoReal_A != PhotoReal_D);
//     EXPECT_TRUE(PhotoReal_A != PhotoReal_E);
//     EXPECT_TRUE(PhotoReal_A != PhotoReal_F);

//     PhotoZeus param_Real_2(mu_def, medium_1, ecuts_1, multiplier_1,
//     hard_component); PhotoBezrukovBugaev param_Real_3(mu_def, medium_1,
//     ecuts_1, multiplier_1, hard_component); PhotoRhode param_Real_4(mu_def,
//     medium_1, ecuts_1, multiplier_1, hard_component); EXPECT_TRUE(PhotoReal_A
//     != param_Real_2); EXPECT_TRUE(PhotoReal_A != param_Real_3);
//     EXPECT_TRUE(PhotoReal_A != param_Real_4);
//     EXPECT_TRUE(param_Real_2 != param_Real_3);
//     EXPECT_TRUE(param_Real_2 != param_Real_4);
//     EXPECT_TRUE(param_Real_3 != param_Real_4);

//     PhotoIntegral Int_PhotoReal_A(PhotoReal_A);
//     PhotoIntegral Int_PhotoReal_B(PhotoReal_B);
//     EXPECT_TRUE(Int_PhotoReal_A != Int_PhotoReal_B);

//     PhotoInterpolant Interpol_PhotoReal_A(PhotoReal_A, InterpolDef);
//     PhotoInterpolant Interpol_PhotoReal_B(PhotoReal_B, InterpolDef);
//     EXPECT_TRUE(Interpol_PhotoReal_A != Interpol_PhotoReal_B);
//     //
//     PhotoAbramowiczLevinLevyMaor97 PhotoQ2_A(mu_def, medium_1, ecuts_1,
//     multiplier_1, shadow_1); PhotoAbramowiczLevinLevyMaor97
//     PhotoQ2_B(tau_def, medium_1, ecuts_1, multiplier_1, shadow_1);
//     PhotoAbramowiczLevinLevyMaor97 PhotoQ2_C(mu_def, medium_2, ecuts_1,
//     multiplier_1, shadow_1); PhotoAbramowiczLevinLevyMaor97 PhotoQ2_D(mu_def,
//     medium_1, ecuts_2, multiplier_1, shadow_1);
//     PhotoAbramowiczLevinLevyMaor97 PhotoQ2_E(mu_def, medium_1, ecuts_1,
//     multiplier_2, shadow_1); EXPECT_TRUE(PhotoQ2_A != PhotoQ2_B);
//     EXPECT_TRUE(PhotoQ2_A != PhotoQ2_C);
//     EXPECT_TRUE(PhotoQ2_A != PhotoQ2_D);
//     EXPECT_TRUE(PhotoQ2_A != PhotoQ2_E);
//     //
//     EXPECT_TRUE(PhotoReal_A != PhotoQ2_A);

//     PhotoAbramowiczLevinLevyMaor91 param_Q2_2(mu_def, medium_1, ecuts_1,
//     multiplier_1, shadow_1); PhotoButkevichMikheyev param_Q2_3(mu_def,
//     medium_1, ecuts_1, multiplier_1, shadow_1); PhotoRenoSarcevicSu
//     param_Q2_4(mu_def, medium_1, ecuts_1, multiplier_1, shadow_1);
//     EXPECT_TRUE(PhotoQ2_A != param_Q2_2);
//     EXPECT_TRUE(PhotoQ2_A != param_Q2_3);
//     EXPECT_TRUE(PhotoQ2_A != param_Q2_4);
//     EXPECT_TRUE(param_Q2_2 != param_Q2_3);
//     EXPECT_TRUE(param_Q2_2 != param_Q2_4);
//     EXPECT_TRUE(param_Q2_3 != param_Q2_4);

//     PhotoIntegral Int_PhotoQ2_A(PhotoQ2_A);
//     PhotoIntegral Int_PhotoQ2_B(PhotoQ2_B);
//     EXPECT_TRUE(Int_PhotoQ2_A != Int_PhotoQ2_B);

//     PhotoQ2Interpolant<PhotoAbramowiczLevinLevyMaor97>
//     PhotoQ2_A_interpol(mu_def, medium_1, ecuts_1, multiplier_1, shadow_1,
//     InterpolDef); PhotoQ2Interpolant<PhotoAbramowiczLevinLevyMaor97>
//     PhotoQ2_B_interpol(tau_def, medium_1, ecuts_1, multiplier_1, shadow_1,
//     InterpolDef); PhotoInterpolant Interpol_PhotoQ2_A(PhotoQ2_A_interpol,
//     InterpolDef); PhotoInterpolant Interpol_PhotoQ2_B(PhotoQ2_B_interpol,
//     InterpolDef); EXPECT_TRUE(Interpol_PhotoQ2_A != Interpol_PhotoQ2_B);
// }

// TEST(Assignment, Copyconstructor)
// {
//     ParticleDef particle_def = MuMinusDef::Get();
//     auto medium = std::make_shared<const Water>();
//     EnergyCutSettings ecuts;
//     double multiplier   = 1.;
//     bool hard_component = true;
//     ShadowButkevichMikheyev shadow;
//     InterpolationDef InterpolDef;

//     PhotoKokoulin PhotoReal_A(particle_def, medium, ecuts, multiplier,
//     hard_component); PhotoKokoulin PhotoReal_B = PhotoReal_A;
//     EXPECT_TRUE(PhotoReal_A == PhotoReal_B);

//     PhotoIntegral Int_PhotoReal_A(PhotoReal_A);
//     PhotoIntegral Int_PhotoReal_B = Int_PhotoReal_A;
//     EXPECT_TRUE(Int_PhotoReal_A == Int_PhotoReal_B);

//     PhotoInterpolant Interpol_PhotoReal_A(PhotoReal_A, InterpolDef);
//     PhotoInterpolant Interpol_PhotoReal_B = Interpol_PhotoReal_A;
//     EXPECT_TRUE(Interpol_PhotoReal_A == Interpol_PhotoReal_B);

//     PhotoAbramowiczLevinLevyMaor97 PhotoQ2_A(particle_def, medium, ecuts,
//     multiplier, shadow); PhotoAbramowiczLevinLevyMaor97 PhotoQ2_B =
//     PhotoQ2_A; EXPECT_TRUE(PhotoQ2_A == PhotoQ2_B);

//     PhotoIntegral Int_PhotoQ2_A(PhotoQ2_A);
//     PhotoIntegral Int_PhotoQ2_B = Int_PhotoQ2_A;
//     EXPECT_TRUE(Int_PhotoQ2_A == Int_PhotoQ2_B);

//     PhotoQ2Interpolant<PhotoAbramowiczLevinLevyMaor97>
//     PhotoQ2_A_interpol(particle_def, medium, ecuts, multiplier, shadow,
//     InterpolDef); PhotoInterpolant Interpol_PhotoQ2_A(PhotoQ2_A_interpol,
//     InterpolDef); PhotoInterpolant Interpol_PhotoQ2_B = Interpol_PhotoQ2_A;
//     EXPECT_TRUE(Interpol_PhotoQ2_A == Interpol_PhotoQ2_B);
// }

// TEST(Assignment, Copyconstructor2)
// {
//     ParticleDef particle_def = MuMinusDef::Get();
//     auto medium = std::make_shared<const Water>();
//     EnergyCutSettings ecuts;
//     double multiplier   = 1.;
//     bool hard_component = true;
//     ShadowButkevichMikheyev shadow;
//     InterpolationDef InterpolDef;

//     PhotoKokoulin PhotoReal_A(particle_def, medium, ecuts, multiplier,
//     hard_component); PhotoKokoulin PhotoReal_B(PhotoReal_A);
//     EXPECT_TRUE(PhotoReal_A == PhotoReal_B);

//     PhotoIntegral Int_PhotoReal_A(PhotoReal_A);
//     PhotoIntegral Int_PhotoReal_B(Int_PhotoReal_A);
//     EXPECT_TRUE(Int_PhotoReal_A == Int_PhotoReal_B);

//     PhotoInterpolant Interpol_PhotoReal_A(PhotoReal_A, InterpolDef);
//     PhotoInterpolant Interpol_PhotoReal_B(Interpol_PhotoReal_A);
//     EXPECT_TRUE(Interpol_PhotoReal_A == Interpol_PhotoReal_B);

//     PhotoAbramowiczLevinLevyMaor97 PhotoQ2_A(particle_def, medium, ecuts,
//     multiplier, shadow); PhotoAbramowiczLevinLevyMaor97 PhotoQ2_B(PhotoQ2_A);
//     EXPECT_TRUE(PhotoQ2_A == PhotoQ2_B);

//     PhotoIntegral Int_PhotoQ2_A(PhotoQ2_A);
//     PhotoIntegral Int_PhotoQ2_B(Int_PhotoQ2_A);
//     EXPECT_TRUE(Int_PhotoQ2_A == Int_PhotoQ2_B);

//     PhotoQ2Interpolant<PhotoAbramowiczLevinLevyMaor97>
//     PhotoQ2_A_interpol(particle_def, medium, ecuts, multiplier, shadow,
//     InterpolDef); PhotoInterpolant Interpol_PhotoQ2_A(PhotoQ2_A_interpol,
//     InterpolDef); PhotoInterpolant Interpol_PhotoQ2_B(Interpol_PhotoQ2_A);
//     EXPECT_TRUE(Interpol_PhotoQ2_A == Interpol_PhotoQ2_B);
// }

// in polymorphism an assignmant and swap operator doesn't make sense

TEST(PhotoRealPhotonAssumption, Test_of_dEdx)
{
    auto in = getTestFiles("Photo_Real_dEdx.txt");

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
        parametrization.erase(0, 5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["hard_component"] = hard_component;

        auto cross = make_photonuclearreal(
            particle_def, *medium, ecuts, false, config);

        dEdx_new = cross->CalculatedEdx(energy) * medium->GetMassDensity();

        EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);
    }
}

TEST(PhotoRealPhotonAssumption, Test_of_dNdx)
{
    auto in = getTestFiles("Photo_Real_dNdx.txt");

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
        parametrization.erase(0, 5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["hard_component"] = hard_component;

        auto cross = make_photonuclearreal(
            particle_def, *medium, ecuts, false, config);

        dNdx_new = cross->CalculatedNdx(energy) * medium->GetMassDensity();

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
    }
}

TEST(PhotoRealPhotonAssumption, Test_of_e)
{
    auto in = getTestFiles("Photo_Real_e.txt");

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
    int failed_ctr = 0;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier
        >> energy >> rnd1 >> rnd2 >> stochastic_loss_stored >> parametrization
        >> hard_component) {
        parametrization.erase(0, 5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["hard_component"] = hard_component;

        auto cross = make_photonuclearreal(
            particle_def, *medium, ecuts, false, config);

        auto dNdx_full = cross->CalculatedNdx(energy);
        auto components = medium->GetComponents();
        double sum = 0;

        for (auto comp : components) {
            double dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());
            sum += dNdx_for_comp;
            if (sum > dNdx_full * rnd2) {
                double rate_new = dNdx_for_comp * rnd1;
                if (ecut == INF and vcut == 1) {
#ifndef NDEBUG
                    EXPECT_DEATH(cross->CalculateStochasticLoss(
                                     comp.GetHash(), energy, rate_new),
                        "");
#endif
                } else {
                    auto v = cross->CalculateStochasticLoss(
                        comp.GetHash(), energy, rate_new);
                    // Results of CalculateStochasticLoss are not identical to
                    // PROPOSAL6 because PROPOSAL6 used random numbers instead
                    // of rates to sample dNdx, which causes different (although
                    // not necessarily less accurate) results in some cases (in
                    // all cases I cross-checked by hand, PROPOSAL7 provided
                    // better results)
                    EXPECT_NEAR(v * energy, stochastic_loss_stored,
                        1E-1 * stochastic_loss_stored);
                    if (std::abs(v * energy - stochastic_loss_stored)
                        > 1e-3 * stochastic_loss_stored)
                        failed_ctr++; // number of cases where this happens
                                      // should be limited

                    // cross check where I am suprised that the accuracy is not
                    // actually better, but since this calculation is not even
                    // used in normal calculations I am fine with it
                    auto rate_calculated
                        = cross->CalculateCumulativeCrosssection(
                            energy, comp.GetHash(), v);
                    EXPECT_NEAR(rate_calculated / dNdx_for_comp, rnd1, 5e-3);
                    break;
                }
            }
        }
    }
    EXPECT_LE(failed_ctr, 100);
}

TEST(PhotoRealPhotonAssumption, Test_of_dEdx_Interpolant)
{
    auto in = getTestFiles("Photo_Real_dEdx.txt");

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
        parametrization.erase(0, 5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["hard_component"] = hard_component;

        auto cross
            = make_photonuclearreal(particle_def, *medium, ecuts, true, config);

        dEdx_new = cross->CalculatedEdx(energy) * medium->GetMassDensity();

        if (hard_component == 1 and energy == 1e5)
            EXPECT_NEAR(dEdx_new, dEdx_stored,
                1e-1 * dEdx_stored); // kink in function (see issue #124)
        else if (vcut * energy == ecut)
            EXPECT_NEAR(
                dEdx_new, dEdx_stored, 1e-1 * dEdx_stored); // kink in function
        else if (parametrization == "Rhode" and energy <= 10000)
            EXPECT_NEAR(dEdx_new, dEdx_stored,
                1e-2 * dEdx_stored); // slight "bumps" in function, hard to
                                     // interpolate
        else
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);
    }
}

TEST(PhotoRealPhotonAssumption, Test_of_dNdx_Interpolant)
{
    auto in = getTestFiles("Photo_Real_dNdx.txt");

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
        parametrization.erase(0, 5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["hard_component"] = hard_component;

        auto cross
            = make_photonuclearreal(particle_def, *medium, ecuts, true, config);

        dNdx_new = cross->CalculatedNdx(energy) * medium->GetMassDensity();
        if (energy * vcut == ecut)
            EXPECT_NEAR(
                dNdx_new, dNdx_stored, 1e-1 * dNdx_stored); // kink in function
        else if (hard_component == 1 and energy >= 1e10)
            EXPECT_NEAR(dNdx_new, dNdx_stored,
                1e-2 * dNdx_stored); // for high E, high_component correction
                                     // produces artefacts due to hard cutoff
                                     // for v=1e-7 in differential cross section
                                     // (see issue #124)
        else if (hard_component == 1 and energy == 1e5)
            EXPECT_NEAR(dNdx_new, dNdx_stored,
                1e-1 * dNdx_stored); // kink in hard_component for E=1e5 (see
                                     // issue #124)
        else
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
    }
}

TEST(PhotoRealPhotonAssumption, Test_of_e_Interpolant)
{
    auto in = getTestFiles("Photo_Real_e.txt");

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
        parametrization.erase(0, 5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

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
        auto dNdx_full = cross->CalculatedNdx(energy);
        auto components = medium->GetComponents();
        double sum = 0;

        for (auto comp : components) {
            double dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());
            sum += dNdx_for_comp;
            if (sum > dNdx_full * rnd2) {
                double rate_new = dNdx_for_comp * rnd1;
                if (ecut == INF and vcut == 1) {
#ifndef NDEBUG
                    EXPECT_DEATH(cross->CalculateStochasticLoss(
                                     comp.GetHash(), energy, rate_new),
                        "");
#endif
                } else {
                    double v = cross->CalculateStochasticLoss(
                        comp.GetHash(), energy, rate_new);
                    // expect correct order of magnitude compared with old
                    // PROPOSAL
                    EXPECT_NEAR(v * energy, stochastic_loss_stored,
                        1.1E-1 * stochastic_loss_stored);

                    // Crosscheck by comparing with integrated value
                    auto rate_calculated
                        = cross_integral->CalculateCumulativeCrosssection(
                            energy, comp.GetHash(), v);
                    EXPECT_NEAR(rate_calculated / dNdx_for_comp, rnd1, 1e-5);

                    break;
                }
            }
        }
    }
}

// Photonuclear parametrization with Q2 Integration

TEST(PhotoQ2Integration, Test_of_dEdx)
{
    auto in = getTestFiles("Photo_Q2_dEdx.txt");

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
        parametrization.erase(0, 5);
        if (parametrization == "ButkevichMikhailov")
            parametrization = "ButkevichMikheyev";
        shadowing.erase(0, 6);
        if (shadowing == "ButkevichMikhailov")
            shadowing = "ButkevichMikheyev";
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["shadow"] = shadowing;

        auto cross
            = make_photonuclearQ2(particle_def, *medium, ecuts, false, config);

        dEdx_new = cross->CalculatedEdx(energy) * medium->GetMassDensity();

        ASSERT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);
    }
}

TEST(PhotoQ2Integration, Test_of_dNdx)
{
    auto in = getTestFiles("Photo_Q2_dNdx.txt");

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
        parametrization.erase(0, 5);
        if (parametrization == "ButkevichMikhailov")
            parametrization = "ButkevichMikheyev";
        shadowing.erase(0, 6);
        if (shadowing == "ButkevichMikhailov")
            shadowing = "ButkevichMikheyev";
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["shadow"] = shadowing;

        auto cross
            = make_photonuclearQ2(particle_def, *medium, ecuts, false, config);

        dNdx_new = cross->CalculatedNdx(energy) * medium->GetMassDensity();

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
    }
}

TEST(PhotoQ2Integration, Test_of_e)
{
    auto in = getTestFiles("Photo_Q2_e.txt");

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
        parametrization.erase(0, 5);
        if (parametrization == "ButkevichMikhailov")
            parametrization = "ButkevichMikheyev";
        shadowing.erase(0, 6);
        if (shadowing == "ButkevichMikhailov")
            shadowing = "ButkevichMikheyev";
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["shadow"] = shadowing;

        auto cross
            = make_photonuclearQ2(particle_def, *medium, ecuts, false, config);

        auto dNdx_full = cross->CalculatedNdx(energy);
        auto components = medium->GetComponents();
        double sum = 0;

        for (auto comp : components) {
            double dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());
            sum += dNdx_for_comp;
            if (sum > dNdx_full * rnd2) {
                double rate_new = dNdx_for_comp * rnd1;
                if (ecut == INF and vcut == 1) {
#ifndef NDEBUG
                    EXPECT_DEATH(cross->CalculateStochasticLoss(
                                     comp.GetHash(), energy, rate_new),
                        "");
#endif
                } else {
                    auto v = cross->CalculateStochasticLoss(
                        comp.GetHash(), energy, rate_new);
                    EXPECT_NEAR(v * energy, stochastic_loss_stored,
                        1E-3 * stochastic_loss_stored);

                    // cross check
                    auto rate_new = cross->CalculateCumulativeCrosssection(
                        energy, comp.GetHash(), v);
                    EXPECT_NEAR(rate_new / dNdx_for_comp, rnd1, 1e-5);
                    break;
                }
            }
        }
    }
}

TEST(PhotoQ2Integration, Test_of_dEdx_Interpolant)
{
    auto in = getTestFiles("Photo_Q2_dEdx.txt");

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
        parametrization.erase(0, 5);
        if (parametrization == "ButkevichMikhailov")
            parametrization = "ButkevichMikheyev";
        shadowing.erase(0, 6);
        if (shadowing == "ButkevichMikhailov")
            shadowing = "ButkevichMikheyev";
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["shadow"] = shadowing;

        auto cross
            = make_photonuclearQ2(particle_def, *medium, ecuts, true, config);

        dEdx_new = cross->CalculatedEdx(energy) * medium->GetMassDensity();
        if (ecut == vcut * energy)
            EXPECT_NEAR(
                dEdx_new, dEdx_stored, 1e-1 * dEdx_stored); // kink in integral
        else if (parametrization == "AbramowiczLevinLevyMaor91"
            and shadowing == "ButkevichMikheyev" and particleName == "MuMinus"
            and energy == 100000000000 and ecut == INF and vcut == 1)
            EXPECT_NEAR(dEdx_new, dEdx_stored,
                1e-2 * dEdx_stored); // bump in integral for these specific
                                     // parameter
        else
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);
    }
}

TEST(PhotoQ2Integration, Test_of_dNdx_Interpolant)
{
    auto in = getTestFiles("Photo_Q2_dNdx.txt");

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
        parametrization.erase(0, 5);
        if (parametrization == "ButkevichMikhailov")
            parametrization = "ButkevichMikheyev";
        shadowing.erase(0, 6);
        if (shadowing == "ButkevichMikhailov")
            shadowing = "ButkevichMikheyev";
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["shadow"] = shadowing;

        auto cross
            = make_photonuclearQ2(particle_def, *medium, ecuts, true, config);

        dNdx_new = cross->CalculatedNdx(energy) * medium->GetMassDensity();
        if (energy * vcut == ecut)
            EXPECT_NEAR(
                dNdx_new, dNdx_stored, 1e-1 * dNdx_stored); // kink in integral
        else
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
    }
}

TEST(PhotoQ2Integration, Test_of_e_Interpolant)
{
    auto in = getTestFiles("Photo_Q2_e.txt");

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
        parametrization.erase(0, 5);
        if (parametrization == "ButkevichMikhailov")
            parametrization = "ButkevichMikheyev";
        shadowing.erase(0, 6);
        if (shadowing == "ButkevichMikhailov")
            shadowing = "ButkevichMikheyev";
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["shadow"] = shadowing;

        auto cross
            = make_photonuclearQ2(particle_def, *medium, ecuts, true, config);
        auto cross_integral
            = make_photonuclearQ2(particle_def, *medium, ecuts, false, config);

        auto dNdx_full = cross->CalculatedNdx(energy);
        auto components = medium->GetComponents();
        double sum = 0;

        for (auto comp : components) {
            double dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());
            sum += dNdx_for_comp;
            if (sum > dNdx_full * rnd2) {
                double rate_new = dNdx_for_comp * rnd1;
                if (ecut == INF and vcut == 1) {
#ifndef NDEBUG
                    EXPECT_DEATH(cross->CalculateStochasticLoss(
                                     comp.GetHash(), energy, rate_new),
                        "");
#endif
                } else {
                    auto v = cross->CalculateStochasticLoss(
                        comp.GetHash(), energy, rate_new);
                    // expect same order of magnitude compared to old PROPOSAL
                    EXPECT_NEAR(energy * v, stochastic_loss_stored,
                        5E-1 * stochastic_loss_stored);

                    // cross check
                    auto rate_rnd
                        = cross_integral->CalculateCumulativeCrosssection(
                            energy, comp.GetHash(), v);
                    if (energy * vcut == ecut or rnd1 < 0.05)
                        EXPECT_NEAR(rate_rnd / dNdx_for_comp, rnd1,
                            1e-2); // kink in integral / TODO: not working too
                                   // well for small rnd
                    else
                        EXPECT_NEAR(rate_rnd / dNdx_for_comp, rnd1, 1e-3);

                    break;
                }
            }
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
