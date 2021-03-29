
#include "gtest/gtest.h"

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

const std::string testfile_dir = "tests/TestFiles/";

// TEST(Comparison, Comparison_equal)
// {
//     ParticleDef particle_def = MuMinusDef::Get();
//     auto medium = std::make_shared<const Water>();
//     EnergyCutSettings ecuts;
//     double multiplier = 1.;
//     bool lpm          = true;

//     EpairProductionRhoIntegral* EpairInt_A =
//         new EpairKelnerKokoulinPetrukhin(particle_def, medium, ecuts, multiplier, lpm);
//     Parametrization* EpairInt_B = new EpairKelnerKokoulinPetrukhin(particle_def, medium, ecuts, multiplier, lpm);
//     EXPECT_TRUE(*EpairInt_A == *EpairInt_B);

//     EpairKelnerKokoulinPetrukhin param_int(particle_def, medium, ecuts, multiplier, lpm);
//     EXPECT_TRUE(param_int == *EpairInt_A);

//     EpairIntegral* Int_A        = new EpairIntegral(param_int);
//     CrossSectionIntegral* Int_B = new EpairIntegral(param_int);
//     EXPECT_TRUE(*Int_A == *Int_B);

//     InterpolationDef InterpolDef;
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin>* EpairInterpol_A =
//         new EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin>(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
//     Parametrization* EpairInterpol_B =
//         new EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin>(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
//     EXPECT_TRUE(*EpairInterpol_A == *EpairInterpol_B);

//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> param_interpol(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
//     EXPECT_TRUE(param_interpol == *EpairInterpol_A);

//     EpairInterpolant* Interpol_A        = new EpairInterpolant(param_interpol, InterpolDef);
//     CrossSectionInterpolant* Interpol_B = new EpairInterpolant(param_interpol, InterpolDef);
//     EXPECT_TRUE(*Interpol_A == *Interpol_B);

//     delete EpairInt_A;
//     delete EpairInt_B;
//     delete Int_A;
//     delete Int_B;
//     delete EpairInterpol_A;
//     delete EpairInterpol_B;
//     delete Interpol_A;
//     delete Interpol_B;
// }

// TEST(Comparison, Comparison_not_equal)
// {
//     ParticleDef mu_def  = MuMinusDef::Get();
//     ParticleDef tau_def = TauMinusDef::Get();
//     auto medium_1 = std::make_shared<const Water>();
//     auto medium_2 = std::make_shared<const Ice>();
//     EnergyCutSettings ecuts_1(500, -1);
//     EnergyCutSettings ecuts_2(-1, 0.05);
//     double multiplier_1 = 1.;
//     double multiplier_2 = 2.;
//     bool lpm_1          = true;
//     bool lpm_2          = false;

//     EpairKelnerKokoulinPetrukhin EpairInt_A(mu_def, medium_1, ecuts_1, multiplier_1, lpm_1);
//     EpairKelnerKokoulinPetrukhin EpairInt_B(tau_def, medium_1, ecuts_1, multiplier_1, lpm_1);
//     EpairKelnerKokoulinPetrukhin EpairInt_C(mu_def, medium_2, ecuts_1, multiplier_1, lpm_1);
//     EpairKelnerKokoulinPetrukhin EpairInt_D(mu_def, medium_1, ecuts_2, multiplier_1, lpm_1);
//     EpairKelnerKokoulinPetrukhin EpairInt_E(mu_def, medium_1, ecuts_1, multiplier_2, lpm_1);
//     EpairKelnerKokoulinPetrukhin EpairInt_F(mu_def, medium_1, ecuts_1, multiplier_1, lpm_2);
//     EXPECT_TRUE(EpairInt_A != EpairInt_B);
//     EXPECT_TRUE(EpairInt_A != EpairInt_C);
//     EXPECT_TRUE(EpairInt_A != EpairInt_D);
//     EXPECT_TRUE(EpairInt_A != EpairInt_E);
//     EXPECT_TRUE(EpairInt_A != EpairInt_F);

//     EpairIntegral Int_A(EpairInt_A);
//     EpairIntegral Int_B(EpairInt_B);
//     EXPECT_TRUE(Int_A != Int_B);

//     InterpolationDef InterpolDef;
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_A(mu_def, medium_1, ecuts_1, multiplier_1, lpm_1, InterpolDef);
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_B(tau_def, medium_1, ecuts_1, multiplier_1, lpm_1, InterpolDef);
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_C(mu_def, medium_2, ecuts_1, multiplier_1, lpm_1, InterpolDef);
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_D(mu_def, medium_1, ecuts_2, multiplier_1, lpm_1, InterpolDef);
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_E(mu_def, medium_1, ecuts_1, multiplier_2, lpm_1, InterpolDef);
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_F(mu_def, medium_1, ecuts_1, multiplier_1, lpm_2, InterpolDef);
//     EXPECT_TRUE(EpairInterpol_A != EpairInterpol_B);
//     EXPECT_TRUE(EpairInterpol_A != EpairInterpol_C);
//     EXPECT_TRUE(EpairInterpol_A != EpairInterpol_D);
//     EXPECT_TRUE(EpairInterpol_A != EpairInterpol_E);
//     EXPECT_TRUE(EpairInterpol_A != EpairInterpol_F);

//     EpairInterpolant Interpol_A(EpairInterpol_A, InterpolDef);
//     EpairInterpolant Interpol_B(EpairInterpol_B, InterpolDef);
//     EXPECT_TRUE(Interpol_A != Interpol_B);
// }

// TEST(Assignment, Copyconstructor)
// {
//     ParticleDef particle_def = MuMinusDef::Get();
//     auto medium = std::make_shared<const Water>();
//     EnergyCutSettings ecuts;
//     double multiplier = 1.;
//     bool lpm          = true;

//     EpairKelnerKokoulinPetrukhin EpairInt_A(particle_def, medium, ecuts, multiplier, lpm);
//     EpairKelnerKokoulinPetrukhin EpairInt_B = EpairInt_A;
//     EXPECT_TRUE(EpairInt_A == EpairInt_B);

//     EpairIntegral Int_A(EpairInt_A);
//     EpairIntegral Int_B = Int_A;
//     EXPECT_TRUE(Int_A == Int_B);

//     InterpolationDef InterpolDef;
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_A(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_B = EpairInterpol_A;
//     EXPECT_TRUE(EpairInterpol_A == EpairInterpol_B);

//     EpairInterpolant Interpol_A(EpairInterpol_A, InterpolDef);
//     EpairInterpolant Interpol_B = Interpol_A;
//     EXPECT_TRUE(Interpol_A == Interpol_B);
// }

// TEST(Assignment, Copyconstructor2)
// {
//     ParticleDef particle_def = MuMinusDef::Get();
//     auto medium = std::make_shared<const Water>();
//     EnergyCutSettings ecuts;
//     double multiplier = 1.;
//     bool lpm          = true;

//     EpairKelnerKokoulinPetrukhin EpairInt_A(particle_def, medium, ecuts, multiplier, lpm);
//     EpairKelnerKokoulinPetrukhin EpairInt_B(EpairInt_A);
//     EXPECT_TRUE(EpairInt_A == EpairInt_B);

//     EpairIntegral Int_A(EpairInt_A);
//     EpairIntegral Int_B(Int_A);
//     EXPECT_TRUE(Int_A == Int_B);

//     InterpolationDef InterpolDef;
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_A(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_B(EpairInterpol_A);
//     EXPECT_TRUE(EpairInterpol_A == EpairInterpol_B);

//     EpairInterpolant Interpol_A(EpairInterpol_A, InterpolDef);
//     EpairInterpolant Interpol_B(Interpol_A);
//     EXPECT_TRUE(Interpol_A == Interpol_B);
// }

// in polymorphism an assignmant and swap operator doesn't make sense

TEST(Epairproduction, Test_of_dEdx)
{
    auto in = getTestFiles("Epair_dEdx.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double dEdx_stored;
    double dEdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> dEdx_stored)
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
        config["lpm"] = lpm;

        auto cross = make_epairproduction(particle_def, *medium, ecuts,
                                          false, config);

        dEdx_new = cross->CalculatedEdx(energy) * medium->GetMassDensity();
        EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-5 * dEdx_stored);
    }
}

TEST(Epairproduction, Test_of_dNdx)
{
    auto in = getTestFiles("Epair_dNdx.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> dNdx_stored)
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
        config["lpm"] = lpm;

        auto cross = make_epairproduction(particle_def, *medium, ecuts, false,
                                          config);

        dNdx_new = cross->CalculatedNdx(energy) * medium->GetMassDensity();

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);
    }
}

TEST(Epairproduction, Test_Stochastic_Loss)
{
    auto in = getTestFiles("Epair_e.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
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

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> rnd1 >> rnd2 >> stochastic_loss_stored)
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
        config["lpm"] = lpm;

        auto cross = make_epairproduction(particle_def, *medium, ecuts, false,
                                          config);

        auto dNdx_full = cross->CalculatedNdx(energy);
        auto components = medium->GetComponents();
        double sum = 0;

        for (auto comp : components)
        {
            double dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());
            sum += dNdx_for_comp;
            if (sum > dNdx_full * rnd2) {
                double rate_new = dNdx_for_comp * rnd1;
                if (ecut == INF and vcut == 1 ) {
                    #ifndef NDEBUG
                    EXPECT_DEATH(cross->CalculateStochasticLoss(comp.GetHash(), energy, rate_new), "");
                    #endif
                } else {
                    auto v =  cross->CalculateStochasticLoss(comp.GetHash(), energy, rate_new);
                    EXPECT_NEAR(v * energy, stochastic_loss_stored, 1E-3 * stochastic_loss_stored);

                    // cross check
                    auto rate_rnd = cross->CalculateCumulativeCrosssection(energy, comp.GetHash(), v);
                    EXPECT_NEAR(rate_rnd/dNdx_for_comp, rnd1, 1e-3);
                    break;
                }
            }
        }
    }
}

TEST(Epairproduction, Test_of_dEdx_Interpolant)
{
    auto in = getTestFiles("Epair_dEdx.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double dEdx_stored;
    double dEdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> dEdx_stored)
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
        config["lpm"] = lpm;

        auto cross = make_epairproduction(particle_def, *medium, ecuts, true,
                                          config);

        dEdx_new = cross->CalculatedEdx(energy) * medium->GetMassDensity();

        if (particleName == "TauMinus" and mediumName == "uranium" and energy == 1e4)
            EXPECT_EQ(dEdx_new, 0.); // lower limit in E for table not precise enough
        else if (vcut * energy == ecut)
            EXPECT_NEAR(dEdx_new, dEdx_stored, 5e-2 * dEdx_stored); // kink in interpolated function
        else if (particleName == "TauMinus" and energy <= 10000)
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-2 * dEdx_stored); // integrand looks bad
        else if (particleName == "TauMinus" and mediumName == "hydrogen" and energy <= 1e5)
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-2 * dEdx_stored); // integrand looks bad
        else
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);
    }
}

TEST(Epairproduction, Test_of_dNdx_Interpolant)
{
    auto in = getTestFiles("Epair_dNdx.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> dNdx_stored)
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
        config["lpm"] = lpm;

        auto cross = make_epairproduction(particle_def, *medium, ecuts, true,
                                          config);

        dNdx_new = cross->CalculatedNdx(energy) * medium->GetMassDensity();
        if (vcut * energy == ecut)
            EXPECT_NEAR(dNdx_new, dNdx_stored, 5e-2 * dNdx_stored);
        else if (particleName == "TauMinus" and mediumName == "hydrogen" and energy <= 1e5)
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-2 * dNdx_stored); // integrand looks bad
        else
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
    }
}

TEST(Epairproduction, Test_of_e_interpol)
{
    auto in = getTestFiles("Epair_e.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
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

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> rnd1 >> rnd2 >> stochastic_loss_stored)
    {        parametrization.erase(0,5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;
        config["lpm"] = lpm;

        auto cross = make_epairproduction(particle_def, *medium, ecuts, true,
                                          config);

        auto dNdx_full = cross->CalculatedNdx(energy);
        auto components = medium->GetComponents();
        double sum = 0;

        for (auto comp : components)
        {
            double dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());
            sum += dNdx_for_comp;
            if (sum > dNdx_full * rnd2) {
                double rate_new = dNdx_for_comp * rnd1;
                if (ecut == INF and vcut == 1 ) {
                    #ifndef NDEBUG
                    EXPECT_DEATH(cross->CalculateStochasticLoss(comp.GetHash(), energy, rate_new), "");
                    #endif
                } else {
                    auto v = cross->CalculateStochasticLoss(comp.GetHash(), energy, rate_new);
                    if (rnd1 < 0.1 or rnd1 > 0.9)
                        EXPECT_NEAR(v * energy, stochastic_loss_stored, 5E-2 * stochastic_loss_stored);
                    else if (energy * vcut == ecut)
                        EXPECT_NEAR(v * energy, stochastic_loss_stored, 1E-2 * stochastic_loss_stored);
                    else if (particleName == "TauMinus" and energy <= 1e5)
                        EXPECT_NEAR(v * energy, stochastic_loss_stored, 1E-2 * stochastic_loss_stored); // integrand problems
                    else if (particleName == "EMinus" and energy >= 1e11)
                        EXPECT_NEAR(v * energy, stochastic_loss_stored, 1E-2 * stochastic_loss_stored);
                    else
                        EXPECT_NEAR(v * energy, stochastic_loss_stored, 1E-3 * stochastic_loss_stored);
                    // cross check
                    auto rate_rnd = cross->CalculateCumulativeCrosssection(energy, comp.GetHash(), v);
                    EXPECT_NEAR(rate_rnd/dNdx_for_comp, rnd1, 1e-3); // this is actually important
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
