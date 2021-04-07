
#include "gtest/gtest.h"

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

const std::string testfile_dir = "tests/TestFiles/";

// TEST(Comparison, Comparison_equal)
// {
//     ParticleDef particle_def = MuMinusDef();
//     auto medium = std::make_shared<const Water>();
//     EnergyCutSettings ecuts;
//     double multiplier = 1.;

//     IonizBetheBlochRossi* Ioniz_A = new IonizBetheBlochRossi(particle_def, medium, ecuts, multiplier);
//     Parametrization* Ioniz_B = new IonizBetheBlochRossi(particle_def, medium, ecuts, multiplier);
//     EXPECT_TRUE(*Ioniz_A == *Ioniz_B);

//     IonizBetheBlochRossi param(particle_def, medium, ecuts, multiplier);
//     EXPECT_TRUE(param == *Ioniz_A);

//     IonizIntegral* Int_A        = new IonizIntegral(param);
//     CrossSectionIntegral* Int_B = new IonizIntegral(param);
//     EXPECT_TRUE(*Int_A == *Int_B);

//     InterpolationDef InterpolDef;
//     IonizInterpolant* Interpol_A        = new IonizInterpolant(param, InterpolDef);
//     CrossSectionInterpolant* Interpol_B = new IonizInterpolant(param, InterpolDef);
//     EXPECT_TRUE(*Interpol_A == *Interpol_B);

//     delete Ioniz_A;
//     delete Ioniz_B;
//     delete Int_A;
//     delete Int_B;
//     delete Interpol_A;
//     delete Interpol_B;
// }

// TEST(Comparison, Comparison_not_equal)
// {
//     ParticleDef mu_def  = MuMinusDef();
//     ParticleDef tau_def = TauMinusDef();
//     auto medium_1 = std::make_shared<const Water>();
//     auto medium_2 = std::make_shared<const Ice>();
//     EnergyCutSettings ecuts_1(500, -1);
//     EnergyCutSettings ecuts_2(-1, 0.05);
//     double multiplier_1 = 1.;
//     double multiplier_2 = 2.;

//     IonizBetheBlochRossi Ioniz_A(mu_def, medium_1, ecuts_1, multiplier_1);
//     IonizBetheBlochRossi Ioniz_B(tau_def, medium_1, ecuts_1, multiplier_1);
//     IonizBetheBlochRossi Ioniz_C(mu_def, medium_2, ecuts_1, multiplier_1);
//     IonizBetheBlochRossi Ioniz_D(mu_def, medium_1, ecuts_2, multiplier_1);
//     IonizBetheBlochRossi Ioniz_E(mu_def, medium_1, ecuts_1, multiplier_2);
//     EXPECT_TRUE(Ioniz_A != Ioniz_B);
//     EXPECT_TRUE(Ioniz_A != Ioniz_C);
//     EXPECT_TRUE(Ioniz_A != Ioniz_D);
//     EXPECT_TRUE(Ioniz_A != Ioniz_E);

//     IonizIntegral Int_A(Ioniz_A);
//     IonizIntegral Int_B(Ioniz_B);
//     EXPECT_TRUE(Int_A != Int_B);

//     InterpolationDef InterpolDef;
//     IonizInterpolant Interpol_A(Ioniz_A, InterpolDef);
//     IonizInterpolant Interpol_B(Ioniz_B, InterpolDef);
//     EXPECT_TRUE(Interpol_A != Interpol_B);
// }

// TEST(Assignment, Copyconstructor)
// {
//     ParticleDef mu_def = MuMinusDef();
//     auto medium = std::make_shared<const Water>();
//     EnergyCutSettings ecuts(500, -1);
//     double multiplier = 1.;

//     IonizBetheBlochRossi Ioniz_A(mu_def, medium, ecuts, multiplier);
//     IonizBetheBlochRossi Ioniz_B = Ioniz_A;

//     IonizIntegral Int_A(Ioniz_A);
//     IonizIntegral Int_B = Int_A;
//     EXPECT_TRUE(Int_A == Int_B);

//     InterpolationDef InterpolDef;
//     IonizInterpolant Interpol_A(Ioniz_A, InterpolDef);
//     IonizInterpolant Interpol_B = Interpol_A;
//     EXPECT_TRUE(Interpol_A == Interpol_B);
// }

// TEST(Assignment, Copyconstructor2)
// {

//     ParticleDef mu_def = MuMinusDef();
//     auto medium = std::make_shared<const Water>();
//     EnergyCutSettings ecuts(500, -1);
//     double multiplier = 1.;

//     IonizBetheBlochRossi Ioniz_A(mu_def, medium, ecuts, multiplier);
//     IonizBetheBlochRossi Ioniz_B(Ioniz_A);

//     IonizIntegral Int_A(Ioniz_A);
//     IonizIntegral Int_B(Int_A);
//     EXPECT_TRUE(Int_A == Int_B);

//     InterpolationDef InterpolDef;
//     IonizInterpolant Interpol_A(Ioniz_A, InterpolDef);
//     IonizInterpolant Interpol_B(Interpol_A);
//     EXPECT_TRUE(Interpol_A == Interpol_B);
// }

// in polymorphism an assignmant and swap operator doesn't make sense

TEST(Ionization, Test_of_dEdx)
{
    auto in = getTestFiles("Ioniz_dEdx.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> dEdx_stored >> parametrization)
    {
        parametrization.erase(0,5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        if (parametrization != "BetheBlochRossi" and particle_def.mass != ME)
            continue;

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_ionization(particle_def, *medium, ecuts, false,
                                     config);

        dEdx_new = cross->CalculatedEdx(energy) * medium->GetMassDensity();
        if (parametrization == "BetheBlochRossi")
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-4 * dEdx_stored); // integration routine changed
        else
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-6 * dEdx_stored);

    }
}

TEST(Ionization, Test_of_dNdx)
{
    auto in = getTestFiles("Ioniz_dNdx.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> dNdx_stored >> parametrization)
    {
        parametrization.erase(0,5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        if (parametrization != "BetheBlochRossi" and particle_def.mass != ME)
            continue;

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_ionization(particle_def, *medium, ecuts, false,
                                     config);

        dNdx_new = cross->CalculatedNdx(energy) * medium->GetMassDensity();
        EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-6 * dNdx_stored);
    }
}

TEST(Ionization, Test_Stochastic_Loss)
{
    auto in = getTestFiles("Ioniz_e.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;
    double stochastic_loss_new;

    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> rnd1 >> rnd2 >> stochastic_loss_stored >> parametrization)
    {
        parametrization.erase(0,5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        if (parametrization != "BetheBlochRossi" and particle_def.mass != ME)
            continue;

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_ionization(particle_def, *medium, ecuts, false,
                                     config);

        auto dNdx_full = cross->CalculatedNdx(energy);
        double sum = 0;

        double dNdx_for_comp = cross->CalculatedNdx(energy, medium->GetHash());
        sum += dNdx_for_comp;
        if (sum > dNdx_full * rnd2) {
            double rate_new = dNdx_for_comp * rnd1;
            if (ecut == INF and vcut == 1 ) {
                #ifndef NDEBUG
                EXPECT_DEATH(cross->CalculateStochasticLoss(medium->GetHash(), energy, rate_new), "");
                #endif
            } else {
                stochastic_loss_new = energy * cross->CalculateStochasticLoss(medium->GetHash(), energy, rate_new);
                EXPECT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);
                break;
            }
        }
    }
}

TEST(Ionization, Test_of_dEdx_Interpolant)
{
    auto in = getTestFiles("Ioniz_dEdx.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> dEdx_stored >> parametrization)
    {
        parametrization.erase(0,5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        if (parametrization != "BetheBlochRossi" and particle_def.mass != ME)
            continue;

        auto cross = make_ionization(particle_def, *medium, ecuts, true,
                                     config);

        dEdx_new = cross->CalculatedEdx(energy) * medium->GetMassDensity();
        if (vcut * energy == ecut)
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-2 * dEdx_stored); // kink in interpolated function in this case
        else
            EXPECT_NEAR(dEdx_new, dEdx_stored, 5e-4 * dEdx_stored);

    }
}

TEST(Ionization, Test_of_dNdx_Interpolant)
{
    auto in = getTestFiles("Ioniz_dNdx.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> dNdx_stored >> parametrization)
    {
        parametrization.erase(0,5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        if (parametrization != "BetheBlochRossi" and particle_def.mass != ME)
            continue;

        auto cross = make_ionization(particle_def, *medium, ecuts, true,
                                     config);

        dNdx_new = cross->CalculatedNdx(energy) * medium->GetMassDensity();
        if (vcut * energy == ecut)
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-1 * dNdx_stored); // kink in interpolated function
        else
            EXPECT_NEAR(dNdx_new, dNdx_stored, 5e-4 * dNdx_stored);
    }
}

TEST(Ionization, Test_of_e_interpol)
{
    auto in = getTestFiles("Ioniz_e.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;
    double stochastic_loss_new;

    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> rnd1 >> rnd2 >> stochastic_loss_stored >> parametrization)
    {
        parametrization.erase(0,5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        if (parametrization != "BetheBlochRossi" and particle_def.mass != ME)
            continue;

        auto cross = make_ionization(particle_def, *medium, ecuts, true,
                                     config);

        auto dNdx_full = cross->CalculatedNdx(energy);
        double sum = 0;

        double dNdx_for_comp = cross->CalculatedNdx(energy, medium->GetHash());
        sum += dNdx_for_comp;
        if (sum > dNdx_full * rnd2) {
            double rate_new = dNdx_for_comp * rnd1;
            if (ecut == INF and vcut == 1 ) {
                #ifndef NDEBUG
                EXPECT_DEATH(cross->CalculateStochasticLoss(medium->GetHash(), energy, rate_new), "");
                #endif
            } else {
                stochastic_loss_new = energy * cross->CalculateStochasticLoss(medium->GetHash(), energy, rate_new);
                EXPECT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-5 * stochastic_loss_stored);
                break;
            }
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
