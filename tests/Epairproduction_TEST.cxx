
#include "gtest/gtest.h"

#include <fstream>
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/EpairIntegral.h"
#include "PROPOSAL/crossection/EpairInterpolant.h"
#include "PROPOSAL/crossection/parametrization/EpairProduction.h"
#include "PROPOSAL/crossection/factories/EpairProductionFactory.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "MuMinus")
    {
        return MuMinusDef::Get();
    } else if (name == "TauMinus")
    {
        return TauMinusDef::Get();
    } else
    {
        return EMinusDef::Get();
    }
}

const std::string testfile_dir = "bin/TestFiles/";

TEST(Comparison, Comparison_equal)
{
    ParticleDef particle_def = MuMinusDef::Get();
    Water medium;
    EnergyCutSettings ecuts;
    double multiplier = 1.;
    bool lpm          = true;

    EpairProductionRhoIntegral* EpairInt_A =
        new EpairKelnerKokoulinPetrukhin(particle_def, medium, ecuts, multiplier, lpm);
    Parametrization* EpairInt_B = new EpairKelnerKokoulinPetrukhin(particle_def, medium, ecuts, multiplier, lpm);
    EXPECT_TRUE(*EpairInt_A == *EpairInt_B);

    EpairKelnerKokoulinPetrukhin param_int(particle_def, medium, ecuts, multiplier, lpm);
    EXPECT_TRUE(param_int == *EpairInt_A);

    EpairIntegral* Int_A        = new EpairIntegral(param_int);
    CrossSectionIntegral* Int_B = new EpairIntegral(param_int);
    EXPECT_TRUE(*Int_A == *Int_B);

    InterpolationDef InterpolDef;
    EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin>* EpairInterpol_A =
        new EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin>(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
    Parametrization* EpairInterpol_B =
        new EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin>(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
    EXPECT_TRUE(*EpairInterpol_A == *EpairInterpol_B);

    EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> param_interpol(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
    EXPECT_TRUE(param_interpol == *EpairInterpol_A);

    EpairInterpolant* Interpol_A        = new EpairInterpolant(param_interpol, InterpolDef);
    CrossSectionInterpolant* Interpol_B = new EpairInterpolant(param_interpol, InterpolDef);
    EXPECT_TRUE(*Interpol_A == *Interpol_B);

    delete EpairInt_A;
    delete EpairInt_B;
    delete Int_A;
    delete Int_B;
    delete EpairInterpol_A;
    delete EpairInterpol_B;
    delete Interpol_A;
    delete Interpol_B;
}

TEST(Comparison, Comparison_not_equal)
{
    ParticleDef mu_def  = MuMinusDef::Get();
    ParticleDef tau_def = TauMinusDef::Get();
    Water medium_1;
    Ice medium_2;
    EnergyCutSettings ecuts_1(500, -1);
    EnergyCutSettings ecuts_2(-1, 0.05);
    double multiplier_1 = 1.;
    double multiplier_2 = 2.;
    bool lpm_1          = true;
    bool lpm_2          = false;

    EpairKelnerKokoulinPetrukhin EpairInt_A(mu_def, medium_1, ecuts_1, multiplier_1, lpm_1);
    EpairKelnerKokoulinPetrukhin EpairInt_B(tau_def, medium_1, ecuts_1, multiplier_1, lpm_1);
    EpairKelnerKokoulinPetrukhin EpairInt_C(mu_def, medium_2, ecuts_1, multiplier_1, lpm_1);
    EpairKelnerKokoulinPetrukhin EpairInt_D(mu_def, medium_1, ecuts_2, multiplier_1, lpm_1);
    EpairKelnerKokoulinPetrukhin EpairInt_E(mu_def, medium_1, ecuts_1, multiplier_2, lpm_1);
    EpairKelnerKokoulinPetrukhin EpairInt_F(mu_def, medium_1, ecuts_1, multiplier_1, lpm_2);
    EXPECT_TRUE(EpairInt_A != EpairInt_B);
    EXPECT_TRUE(EpairInt_A != EpairInt_C);
    EXPECT_TRUE(EpairInt_A != EpairInt_D);
    EXPECT_TRUE(EpairInt_A != EpairInt_E);
    EXPECT_TRUE(EpairInt_A != EpairInt_F);

    EpairIntegral Int_A(EpairInt_A);
    EpairIntegral Int_B(EpairInt_B);
    EXPECT_TRUE(Int_A != Int_B);

    InterpolationDef InterpolDef;
    EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_A(mu_def, medium_1, ecuts_1, multiplier_1, lpm_1, InterpolDef);
    EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_B(tau_def, medium_1, ecuts_1, multiplier_1, lpm_1, InterpolDef);
    EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_C(mu_def, medium_2, ecuts_1, multiplier_1, lpm_1, InterpolDef);
    EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_D(mu_def, medium_1, ecuts_2, multiplier_1, lpm_1, InterpolDef);
    EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_E(mu_def, medium_1, ecuts_1, multiplier_2, lpm_1, InterpolDef);
    EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_F(mu_def, medium_1, ecuts_1, multiplier_1, lpm_2, InterpolDef);
    EXPECT_TRUE(EpairInterpol_A != EpairInterpol_B);
    EXPECT_TRUE(EpairInterpol_A != EpairInterpol_C);
    EXPECT_TRUE(EpairInterpol_A != EpairInterpol_D);
    EXPECT_TRUE(EpairInterpol_A != EpairInterpol_E);
    EXPECT_TRUE(EpairInterpol_A != EpairInterpol_F);

    EpairInterpolant Interpol_A(EpairInterpol_A, InterpolDef);
    EpairInterpolant Interpol_B(EpairInterpol_B, InterpolDef);
    EXPECT_TRUE(Interpol_A != Interpol_B);
}

TEST(Assignment, Copyconstructor)
{
    ParticleDef particle_def = MuMinusDef::Get();
    Water medium;
    EnergyCutSettings ecuts;
    double multiplier = 1.;
    bool lpm          = true;

    EpairKelnerKokoulinPetrukhin EpairInt_A(particle_def, medium, ecuts, multiplier, lpm);
    EpairKelnerKokoulinPetrukhin EpairInt_B = EpairInt_A;
    EXPECT_TRUE(EpairInt_A == EpairInt_B);

    EpairIntegral Int_A(EpairInt_A);
    EpairIntegral Int_B = Int_A;
    EXPECT_TRUE(Int_A == Int_B);

    InterpolationDef InterpolDef;
    EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_A(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
    EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_B = EpairInterpol_A;
    EXPECT_TRUE(EpairInterpol_A == EpairInterpol_B);

    EpairInterpolant Interpol_A(EpairInterpol_A, InterpolDef);
    EpairInterpolant Interpol_B = Interpol_A;
    EXPECT_TRUE(Interpol_A == Interpol_B);
}

TEST(Assignment, Copyconstructor2)
{
    ParticleDef particle_def = MuMinusDef::Get();
    Water medium;
    EnergyCutSettings ecuts;
    double multiplier = 1.;
    bool lpm          = true;

    EpairKelnerKokoulinPetrukhin EpairInt_A(particle_def, medium, ecuts, multiplier, lpm);
    EpairKelnerKokoulinPetrukhin EpairInt_B(EpairInt_A);
    EXPECT_TRUE(EpairInt_A == EpairInt_B);

    EpairIntegral Int_A(EpairInt_A);
    EpairIntegral Int_B(Int_A);
    EXPECT_TRUE(Int_A == Int_B);

    InterpolationDef InterpolDef;
    EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_A(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
    EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_B(EpairInterpol_A);
    EXPECT_TRUE(EpairInterpol_A == EpairInterpol_B);

    EpairInterpolant Interpol_A(EpairInterpol_A, InterpolDef);
    EpairInterpolant Interpol_B(Interpol_A);
    EXPECT_TRUE(Interpol_A == Interpol_B);
}

// in polymorphism an assignmant and swap operator doesn't make sense

TEST(Epairproduction, Test_of_dEdx)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Epair_dEdx.txt";
    in.open(filename.c_str());

    if (!in.good())
    {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double dEdx_stored;
    double dEdx_new;

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> dEdx_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        EpairProductionFactory::Definition epair_def;
        epair_def.multiplier      = multiplier;
        epair_def.parametrization = EpairProductionFactory::Get().GetEnumFromString(parametrization);
        epair_def.lpm_effect  = lpm;

        CrossSection* epair = EpairProductionFactory::Get().CreateEpairProduction(particle_def, *medium, ecuts, epair_def);
        dEdx_new = epair->CalculatedEdx(energy);

        ASSERT_NEAR(dEdx_new, dEdx_stored, 1e-10 * dEdx_stored);

        delete medium;
        delete epair;
    }
}

TEST(Epairproduction, Test_of_dNdx)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Epair_dNdx.txt";
    in.open(filename.c_str());

    if (!in.good())
    {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> dNdx_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        EpairProductionFactory::Definition epair_def;
        epair_def.multiplier      = multiplier;
        epair_def.parametrization = EpairProductionFactory::Get().GetEnumFromString(parametrization);
        epair_def.lpm_effect  = lpm;

        CrossSection* epair = EpairProductionFactory::Get().CreateEpairProduction(particle_def, *medium, ecuts, epair_def);
        dNdx_new = epair->CalculatedNdx(energy);

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);

        delete medium;
        delete epair;
    }
}

TEST(Epairproduction, Test_of_dNdx_rnd)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Epair_dNdx_rnd.txt";
    in.open(filename.c_str());

    if (!in.good())
    {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double rnd;
    double dNdx_rnd_stored;
    double dNdx_rnd_new;

    RandomGenerator::Get().SetSeed(0);

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> rnd >> dNdx_rnd_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        EpairProductionFactory::Definition epair_def;
        epair_def.multiplier      = multiplier;
        epair_def.parametrization = EpairProductionFactory::Get().GetEnumFromString(parametrization);
        epair_def.lpm_effect  = lpm;

        CrossSection* epair = EpairProductionFactory::Get().CreateEpairProduction(particle_def, *medium, ecuts, epair_def);

        dNdx_rnd_new = epair->CalculatedNdx(energy, rnd);

        ASSERT_NEAR(dNdx_rnd_new, dNdx_rnd_stored, 1E-10 * dNdx_rnd_stored);

        delete medium;
        delete epair;
    }
}

TEST(Epairproduction, Test_Stochastic_Loss)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Epair_e.txt";
    in.open(filename.c_str());

    if (!in.good())
    {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double rnd1, rnd2;
    double stochastic_loss_stored;
    double stochastic_loss_new;

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(0);

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> rnd1 >> rnd2 >>
            stochastic_loss_stored;


        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        EpairProductionFactory::Definition epair_def;
        epair_def.multiplier      = multiplier;
        epair_def.parametrization = EpairProductionFactory::Get().GetEnumFromString(parametrization);
        epair_def.lpm_effect  = lpm;

        CrossSection* epair = EpairProductionFactory::Get().CreateEpairProduction(particle_def, *medium, ecuts, epair_def);

        stochastic_loss_new = epair->CalculateStochasticLoss(energy, rnd1, rnd2);

        ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);

        delete medium;
        delete epair;
    }
}

TEST(Epairproduction, Test_of_dEdx_Interpolant)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Epair_dEdx_interpol.txt";
    in.open(filename.c_str());

    if (!in.good())
    {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double dEdx_stored;
    double dEdx_new;

    InterpolationDef InterpolDef;

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> dEdx_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        EpairProductionFactory::Definition epair_def;
        epair_def.multiplier      = multiplier;
        epair_def.parametrization = EpairProductionFactory::Get().GetEnumFromString(parametrization);
        epair_def.lpm_effect  = lpm;

        CrossSection* epair = EpairProductionFactory::Get().CreateEpairProduction(particle_def, *medium, ecuts, epair_def, InterpolDef);

        dEdx_new = epair->CalculatedEdx(energy);

        ASSERT_NEAR(dEdx_new, dEdx_stored, 1e-10 * dEdx_stored);

        delete medium;
        delete epair;
    }
}

TEST(Epairproduction, Test_of_dNdx_Interpolant)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Epair_dNdx_interpol.txt";
    in.open(filename.c_str());

    if (!in.good())
    {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    InterpolationDef InterpolDef;

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> dNdx_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        EpairProductionFactory::Definition epair_def;
        epair_def.multiplier      = multiplier;
        epair_def.parametrization = EpairProductionFactory::Get().GetEnumFromString(parametrization);
        epair_def.lpm_effect  = lpm;

        CrossSection* epair = EpairProductionFactory::Get().CreateEpairProduction(particle_def, *medium, ecuts, epair_def, InterpolDef);

        dNdx_new = epair->CalculatedNdx(energy);

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);

        delete medium;
        delete epair;
    }
}

TEST(Epairproduction, Test_of_dNdxrnd_interpol)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Epair_dNdx_rnd_interpol.txt";
    in.open(filename.c_str());

    if (!in.good())
    {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double rnd;
    double dNdx_rnd_stored;
    double dNdx_rnd_new;

    InterpolationDef InterpolDef;

    RandomGenerator::Get().SetSeed(0);

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> rnd >> dNdx_rnd_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        EpairProductionFactory::Definition epair_def;
        epair_def.multiplier      = multiplier;
        epair_def.parametrization = EpairProductionFactory::Get().GetEnumFromString(parametrization);
        epair_def.lpm_effect  = lpm;

        CrossSection* epair = EpairProductionFactory::Get().CreateEpairProduction(particle_def, *medium, ecuts, epair_def, InterpolDef);

        dNdx_rnd_new = epair->CalculatedNdx(energy, rnd);

        ASSERT_NEAR(dNdx_rnd_new, dNdx_rnd_stored, 1E-10 * dNdx_rnd_stored);

        delete medium;
        delete epair;
    }
}

TEST(Epairproduction, Test_of_e_interpol)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Epair_e_interpol.txt";
    in.open(filename.c_str());

    if (!in.good())
    {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double rnd1, rnd2;
    double stochastic_loss_stored;
    double stochastic_loss_new;

    InterpolationDef InterpolDef;

    RandomGenerator::Get().SetSeed(0);

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> rnd1 >> rnd2 >>
            stochastic_loss_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        EpairProductionFactory::Definition epair_def;
        epair_def.multiplier      = multiplier;
        epair_def.parametrization = EpairProductionFactory::Get().GetEnumFromString(parametrization);
        epair_def.lpm_effect  = lpm;

        CrossSection* epair = EpairProductionFactory::Get().CreateEpairProduction(particle_def, *medium, ecuts, epair_def, InterpolDef);

        stochastic_loss_new = epair->CalculateStochasticLoss(energy, rnd1, rnd2);

        ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);

        delete medium;
        delete epair;
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
