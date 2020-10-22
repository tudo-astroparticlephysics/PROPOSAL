
#include "gtest/gtest.h"

#include <fstream>
#include "PROPOSAL/crosssection/parametrization/Annihilation.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"

using namespace PROPOSAL;

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "MuMinus")
    {
        return MuMinusDef();
    } else if (name == "TauMinus")
    {
        return TauMinusDef();
    } else if (name == "EMinus")
    {
        return EMinusDef();
    } else if (name == "MuPlus")
    {
        return MuPlusDef();
    } else if (name == "TauPlus")
    {
        return TauPlusDef();
    } else if (name == "EPlus")
    {
        return EPlusDef();
    }
    else{
        return MuMinusDef();
    }
}

// TEST(Comparison, Comparison_equal_particle)
// {
// ParticleDef particle_def = EPlusDef();
// auto medium = std::make_shared<const Water>();
// double multiplier   = 1.;

// Annihilation* Anni_A = new AnnihilationHeitler(particle_def, medium, multiplier);
// Annihilation* Anni_B = new AnnihilationHeitler(particle_def, medium, multiplier);
// EXPECT_TRUE(*Anni_A == *Anni_B);

// AnnihilationHeitler param_int(particle_def, medium, multiplier);
// EXPECT_TRUE(param_int == *Anni_A);

// AnnihilationIntegral* Int_A        = new AnnihilationIntegral(param_int);
// CrossSectionIntegral* Int_B = new AnnihilationIntegral(param_int);
// EXPECT_TRUE(*Int_A == *Int_B);

// InterpolationDef InterpolDef;

// AnnihilationInterpolant* Interpol_A        = new AnnihilationInterpolant(param_int, InterpolDef);
// CrossSectionInterpolant* Interpol_B = new AnnihilationInterpolant(param_int, InterpolDef);
// EXPECT_TRUE(*Interpol_A == *Interpol_B);

// delete Anni_A;
// delete Anni_B;
// delete Int_A;
// delete Int_B;
// delete Interpol_A;
// delete Interpol_B;
// }

// TEST(Comparison, Comparison_not_equal)
// {
// ParticleDef mu_def  = MuMinusDef();
// ParticleDef e_minus_def = EMinusDef();
// ParticleDef e_plus_def  = EPlusDef();
// auto medium_1 = std::make_shared<const Water>();
// auto medium_2 = std::make_shared<const Ice>();
// double multiplier_1 = 1.;
// double multiplier_2 = 2.;

// AnnihilationHeitler Anni_A(mu_def, medium_1, multiplier_1);
// AnnihilationHeitler Anni_B(e_minus_def, medium_1, multiplier_1);
// AnnihilationHeitler Anni_C(mu_def, medium_2, multiplier_1);
// AnnihilationHeitler Anni_D(e_plus_def, medium_1, multiplier_1);
// AnnihilationHeitler Anni_E(mu_def, medium_1, multiplier_2);
// EXPECT_TRUE(Anni_A != Anni_B);
// EXPECT_TRUE(Anni_A != Anni_C);
// EXPECT_TRUE(Anni_A != Anni_D);
// EXPECT_TRUE(Anni_A != Anni_E);

// AnnihilationIntegral Int_A(Anni_A);
// AnnihilationIntegral Int_B(Anni_B);
// EXPECT_TRUE(Int_A != Int_B);

// InterpolationDef InterpolDef;
// AnnihilationIntegral AnniIntegral_A(Anni_A);
// AnnihilationIntegral AnniIntegral_B(Anni_B);
// AnnihilationIntegral AnniIntegral_C(Anni_C);
// AnnihilationIntegral AnniIntegral_D(Anni_D);
// AnnihilationIntegral AnniIntegral_E(Anni_E);
// EXPECT_TRUE(AnniIntegral_A != AnniIntegral_B);
// EXPECT_TRUE(AnniIntegral_A != AnniIntegral_C);
// EXPECT_TRUE(AnniIntegral_A != AnniIntegral_D);
// EXPECT_TRUE(AnniIntegral_A != AnniIntegral_E);

// AnnihilationIntegral Integral_A(AnniIntegral_A);
// AnnihilationIntegral Integral_B(AnniIntegral_B);
// EXPECT_TRUE(Integral_A != Integral_B);
// }

// TEST(Assignment, Copyconstructor)
// {
// ParticleDef particle_def = EPlusDef();
// auto medium = std::make_shared<const Water>();
// double multiplier = 1.;

// AnnihilationHeitler Anni_A(particle_def, medium, multiplier);
// AnnihilationHeitler Anni_B = Anni_A;
// EXPECT_TRUE(Anni_A == Anni_B);

// AnnihilationIntegral Int_A(Anni_A);
// AnnihilationIntegral Int_B = Int_A;
// EXPECT_TRUE(Int_A == Int_B);

// InterpolationDef InterpolDef;
// AnnihilationInterpolant AnniInterpol_A(Anni_A, InterpolDef);
// AnnihilationInterpolant AnniInterpol_B = AnniInterpol_A;
// EXPECT_TRUE(AnniInterpol_A == AnniInterpol_B);
// }




TEST(Annihilation, Test_of_dNdx) {

    std::string filename = "bin/TestFiles/Anni_dNdx.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";


    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double multiplier;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    while (in.good())
    {
        in >> particleName >> mediumName >> multiplier >> energy >> parametrization >> dNdx_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        auto cross = crosssection::make_annihilation(particle_def,
            *medium,
            false,
            parametrization);

        dNdx_new = cross->CalculatedNdx(energy);

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);
    }
}


TEST(Annihilation, Test_Stochastic_Loss)
{
    std::string filename = "bin/TestFiles/Anni_e.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double multiplier;
    double energy;
    std::string parametrization;
    double rnd;
    double rate;
    double stochastic_loss_stored;
    double stochastic_loss_new;

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(0);

    while (in.good())
    {
        in >> particleName >> mediumName >> multiplier >> energy >> parametrization >> rnd >>
        stochastic_loss_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        auto cross = crosssection::make_annihilation(particle_def,
            *medium,
            false,
            parametrization);

        auto components = medium->GetComponents();
        for (auto comp : components)
        {
            auto comp_ptr = std::make_shared<const Component>(comp);
            // first calculate the complete rate, then sample the loss to a rate
            rate = cross->CalculatedNdx(energy, comp_ptr);
            stochastic_loss_new = cross->CalculateStochasticLoss(
                comp_ptr, energy, rnd*rate);

            ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);
        }
    }
}


TEST(Annihilation, Test_of_dNdx_Interpolant)
{
    std::string filename = "bin/TestFiles/Anni_dNdx_interpol.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double multiplier;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    while (in.good())
    {
        in >> particleName >> mediumName >> multiplier >> energy >> parametrization >> dNdx_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        auto cross = crosssection::make_annihilation(particle_def,
            *medium,
            true,
            parametrization);

        dNdx_new = cross->CalculatedNdx(energy);

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);
    }
}

TEST(Annihilation, Test_of_e_interpol)
{
    std::string filename = "bin/TestFiles/Anni_e_interpol.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double multiplier;
    double energy;
    std::string parametrization;
    double rnd;
    double rate;
    double stochastic_loss_stored;
    double stochastic_loss_new;

    RandomGenerator::Get().SetSeed(0);

    while (in.good())
    {
        in >> particleName >> mediumName >> multiplier >>  energy >> parametrization >> rnd >>
        stochastic_loss_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        auto cross = crosssection::make_annihilation(particle_def,
            *medium,
            true,
            parametrization);

        auto components = medium->GetComponents();
        for (auto comp : components)
        {
            auto comp_ptr = std::make_shared<const Component>(comp);
            // first calculate the complete rate, then sample the loss to a rate
            rate = cross->CalculatedNdx(energy, comp_ptr);
            stochastic_loss_new = cross->CalculateStochasticLoss(
                comp_ptr, energy, rnd*rate);

            ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
