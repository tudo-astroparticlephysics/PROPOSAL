
#include "gtest/gtest.h"

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/Factories/WeakInteractionFactory.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSALTestUtilities/TestFilesHandling.h"

using namespace PROPOSAL;

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "MuMinus") {
        return MuMinusDef();
    } else if (name == "TauMinus") {
        return TauMinusDef();
    } else if (name == "EMinus") {
        return EMinusDef();
    } else if (name == "MuPlus") {
        return MuPlusDef();
    } else if (name == "TauPlus") {
        return TauPlusDef();
    } else if (name == "EPlus") {
        return EPlusDef();
    } else {
        return MuMinusDef();
    }
}

// TEST(Comparison, Comparison_equal_particle)
// {
// ParticleDef particle_def = MuMinusDef::Get(); //particle
// auto medium = std::make_shared<const Water>();
// EnergyCutSettings ecuts;
// double multiplier   = 1.;

// WeakInteraction* Weak_A = new WeakCooperSarkarMertsch(particle_def, medium,
// multiplier); Parametrization* Weak_B = new
// WeakCooperSarkarMertsch(particle_def, medium, multiplier);
// EXPECT_TRUE(*Weak_A == *Weak_B);

// WeakCooperSarkarMertsch param_int(particle_def, medium, multiplier);
// EXPECT_TRUE(param_int == *Weak_A);

// WeakIntegral* Int_A        = new WeakIntegral(param_int);
// CrossSectionIntegral* Int_B = new WeakIntegral(param_int);
// EXPECT_TRUE(*Int_A == *Int_B);

// InterpolationDef InterpolDef;

// WeakInterpolant* Interpol_A        = new WeakInterpolant(param_int,
// InterpolDef); CrossSectionInterpolant* Interpol_B = new
// WeakInterpolant(param_int, InterpolDef); EXPECT_TRUE(*Interpol_A ==
// *Interpol_B);

// delete Weak_A;
// delete Weak_B;
// delete Int_A;
// delete Int_B;
// delete Interpol_A;
// delete Interpol_B;
// }

// TEST(Comparison, Comparison_equal_antiparticle)
// {
//     ParticleDef particle_def = MuPlusDef::Get(); //antiparticle
//     auto medium = std::make_shared<const Water>();
//     double multiplier   = 1.;

//     WeakInteraction* Weak_A = new WeakCooperSarkarMertsch(particle_def,
//     medium, multiplier); Parametrization* Weak_B = new
//     WeakCooperSarkarMertsch(particle_def, medium, multiplier);
//     EXPECT_TRUE(*Weak_B == *Weak_B);

//     WeakCooperSarkarMertsch param_int(particle_def, medium, multiplier);
//     EXPECT_TRUE(param_int == *Weak_A);

//     WeakIntegral* Int_A        = new WeakIntegral(param_int);
//     CrossSectionIntegral* Int_B = new WeakIntegral(param_int);
//     EXPECT_TRUE(*Int_A == *Int_B);

//     InterpolationDef InterpolDef;

//     WeakInterpolant* Interpol_A        = new WeakInterpolant(param_int,
//     InterpolDef); CrossSectionInterpolant* Interpol_B = new
//     WeakInterpolant(param_int, InterpolDef); EXPECT_TRUE(*Interpol_A ==
//     *Interpol_B);

//     delete Weak_A;
//     delete Weak_B;
//     delete Int_A;
//     delete Int_B;
//     delete Interpol_A;
//     delete Interpol_B;
// }

// TEST(Comparison, Comparison_not_equal)
// {
// ParticleDef mu_def  = MuMinusDef::Get();
// ParticleDef tau_def = TauMinusDef::Get();
// ParticleDef mu_plus_def  = MuPlusDef::Get();
// auto medium_1 = std::make_shared<const Water>();
// auto medium_2 = std::make_shared<const Ice>();
// double multiplier_1 = 1.;
// double multiplier_2 = 2.;

// WeakCooperSarkarMertsch Weak_A(mu_def, medium_1, multiplier_1);
// WeakCooperSarkarMertsch Weak_B(tau_def, medium_1, multiplier_1);
// WeakCooperSarkarMertsch Weak_C(mu_def, medium_2, multiplier_1);
// WeakCooperSarkarMertsch Weak_D(mu_plus_def, medium_1, multiplier_1);
// WeakCooperSarkarMertsch Weak_E(mu_def, medium_1, multiplier_2);
// EXPECT_TRUE(Weak_A != Weak_B);
// EXPECT_TRUE(Weak_A != Weak_C);
// EXPECT_TRUE(Weak_A != Weak_D);
// EXPECT_TRUE(Weak_A != Weak_E);

// WeakIntegral Int_A(Weak_A);
// WeakIntegral Int_B(Weak_B);
// EXPECT_TRUE(Int_A != Int_B);

// InterpolationDef InterpolDef;
// WeakIntegral WeakIntegral_A(Weak_A);
// WeakIntegral WeakIntegral_B(Weak_B);
// WeakIntegral WeakIntegral_C(Weak_C);
// WeakIntegral WeakIntegral_D(Weak_D);
// WeakIntegral WeakIntegral_E(Weak_E);
// EXPECT_TRUE(WeakIntegral_A != WeakIntegral_B);
// EXPECT_TRUE(WeakIntegral_A != WeakIntegral_C);
// EXPECT_TRUE(WeakIntegral_A != WeakIntegral_D);
// EXPECT_TRUE(WeakIntegral_A != WeakIntegral_E);

// WeakIntegral Integral_A(WeakIntegral_A);
// WeakIntegral Integral_B(WeakIntegral_B);
// EXPECT_TRUE(Integral_A != Integral_B);
// }

// TEST(Assignment, Copyconstructor)
// {
// ParticleDef particle_def = MuMinusDef::Get();
// auto medium = std::make_shared<const Water>();
// double multiplier = 1.;

// WeakCooperSarkarMertsch Weak_A(particle_def, medium, multiplier);
// WeakCooperSarkarMertsch Weak_B = Weak_A;
// EXPECT_TRUE(Weak_A == Weak_B);

// WeakIntegral Int_A(Weak_A);
// WeakIntegral Int_B = Int_A;
// EXPECT_TRUE(Int_A == Int_B);

// InterpolationDef InterpolDef;
// WeakInterpolant WeakInterpol_A(Weak_A, InterpolDef);
// WeakInterpolant WeakInterpol_B = WeakInterpol_A;
// EXPECT_TRUE(WeakInterpol_A == WeakInterpol_B);
// }

TEST(WeakInteraction, Test_of_dNdx)
{
    auto in = getTestFiles("Weak_dNdx.txt");

    std::string particleName;
    std::string mediumName;
    double multiplier;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> multiplier >> energy
        >> parametrization >> dNdx_stored) {
        parametrization.erase(0, 4);
        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_weakinteraction(particle_def, *medium, false, config);

        dNdx_new = cross->CalculatedNdx(energy) * medium->GetMassDensity();

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);
    }
}

TEST(WeakInteraction, Test_Stochastic_Loss)
{
    auto in = getTestFiles("Weak_e.txt");

    std::string particleName;
    std::string mediumName;
    double multiplier;
    double energy;
    std::string parametrization;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;
    double stochastic_loss_new;

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> multiplier >> energy
        >> parametrization >> rnd1 >> rnd2 >> stochastic_loss_stored) {
        parametrization.erase(0, 4);
        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_weakinteraction(particle_def, *medium, false, config);

        auto dNdx_full = cross->CalculatedNdx(energy);
        auto components = medium->GetComponents();
        double sum = 0;

        for (auto comp : components) {
            double dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());
            sum += dNdx_for_comp;
            if (sum > dNdx_full * rnd2) {
                double rate_new = dNdx_for_comp * rnd1;
                stochastic_loss_new = energy
                    * cross->CalculateStochasticLoss(
                        comp.GetHash(), energy, rate_new);
                EXPECT_NEAR(stochastic_loss_new, stochastic_loss_stored,
                    1E-3 * stochastic_loss_stored);
                break;
            }
        }
    }
}

TEST(WeakInteraction, Test_of_dNdx_Interpolant)
{
    auto in = getTestFiles("Weak_dNdx.txt");

    std::string particleName;
    std::string mediumName;
    double multiplier;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> multiplier >> energy
        >> parametrization >> dNdx_stored) {
        parametrization.erase(0, 4);
        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_weakinteraction(particle_def, *medium, true, config);

        dNdx_new = cross->CalculatedNdx(energy) * medium->GetMassDensity();

        EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
    }
}

TEST(WeakInteraction, Test_of_e_interpol)
{
    auto in = getTestFiles("Weak_e.txt");

    std::string particleName;
    std::string mediumName;
    double multiplier;
    double energy;
    std::string parametrization;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;

    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> multiplier >> energy
        >> parametrization >> rnd1 >> rnd2 >> stochastic_loss_stored) {
        parametrization.erase(0, 4);
        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_weakinteraction(particle_def, *medium, true, config);

        auto dNdx_full = cross->CalculatedNdx(energy);
        auto components = medium->GetComponents();
        double sum = 0;
        for (auto comp : components) {
            double dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());
            sum += dNdx_for_comp;
            if (sum >= dNdx_full * rnd2) {
                double rate_new = dNdx_for_comp * rnd1;
                auto v = cross->CalculateStochasticLoss(
                    comp.GetHash(), energy, rate_new);
                EXPECT_NEAR(energy * v, stochastic_loss_stored,
                    5e-2 * stochastic_loss_stored);

                // cross check
                auto rate_rnd = cross->CalculateCumulativeCrosssection(
                    energy, comp.GetHash(), v);
                EXPECT_NEAR(rate_rnd / dNdx_for_comp, rnd1, 1e-4);
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
