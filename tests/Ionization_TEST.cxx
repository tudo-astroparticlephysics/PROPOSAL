
#include "gtest/gtest.h"

#include <fstream>
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/IonizIntegral.h"
#include "PROPOSAL/crossection/IonizInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"
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

    Ionization* Ioniz_A      = new Ionization(particle_def, medium, ecuts, multiplier);
    Parametrization* Ioniz_B = new Ionization(particle_def, medium, ecuts, multiplier);
    EXPECT_TRUE(*Ioniz_A == *Ioniz_B);

    Ionization param(particle_def, medium, ecuts, multiplier);
    EXPECT_TRUE(param == *Ioniz_A);

    IonizIntegral* Int_A        = new IonizIntegral(param);
    CrossSectionIntegral* Int_B = new IonizIntegral(param);
    EXPECT_TRUE(*Int_A == *Int_B);

    InterpolationDef InterpolDef;
    IonizInterpolant* Interpol_A        = new IonizInterpolant(param, InterpolDef);
    CrossSectionInterpolant* Interpol_B = new IonizInterpolant(param, InterpolDef);
    EXPECT_TRUE(*Interpol_A == *Interpol_B);

    delete Ioniz_A;
    delete Ioniz_B;
    delete Int_A;
    delete Int_B;
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

    Ionization Ioniz_A(mu_def, medium_1, ecuts_1, multiplier_1);
    Ionization Ioniz_B(tau_def, medium_1, ecuts_1, multiplier_1);
    Ionization Ioniz_C(mu_def, medium_2, ecuts_1, multiplier_1);
    Ionization Ioniz_D(mu_def, medium_1, ecuts_2, multiplier_1);
    Ionization Ioniz_E(mu_def, medium_1, ecuts_1, multiplier_2);
    EXPECT_TRUE(Ioniz_A != Ioniz_B);
    EXPECT_TRUE(Ioniz_A != Ioniz_C);
    EXPECT_TRUE(Ioniz_A != Ioniz_D);
    EXPECT_TRUE(Ioniz_A != Ioniz_E);

    IonizIntegral Int_A(Ioniz_A);
    IonizIntegral Int_B(Ioniz_B);
    EXPECT_TRUE(Int_A != Int_B);

    InterpolationDef InterpolDef;
    IonizInterpolant Interpol_A(Ioniz_A, InterpolDef);
    IonizInterpolant Interpol_B(Ioniz_B, InterpolDef);
    EXPECT_TRUE(Interpol_A != Interpol_B);
}

TEST(Assignment, Copyconstructor)
{
    ParticleDef mu_def = MuMinusDef::Get();
    Water medium;
    EnergyCutSettings ecuts(500, -1);
    double multiplier = 1.;

    Ionization Ioniz_A(mu_def, medium, ecuts, multiplier);
    Ionization Ioniz_B = Ioniz_A;

    IonizIntegral Int_A(Ioniz_A);
    IonizIntegral Int_B = Int_A;
    EXPECT_TRUE(Int_A == Int_B);

    InterpolationDef InterpolDef;
    IonizInterpolant Interpol_A(Ioniz_A, InterpolDef);
    IonizInterpolant Interpol_B = Interpol_A;
    EXPECT_TRUE(Interpol_A == Interpol_B);
}

TEST(Assignment, Copyconstructor2)
{

    ParticleDef mu_def = MuMinusDef::Get();
    Water medium;
    EnergyCutSettings ecuts(500, -1);
    double multiplier = 1.;

    Ionization Ioniz_A(mu_def, medium, ecuts, multiplier);
    Ionization Ioniz_B(Ioniz_A);

    IonizIntegral Int_A(Ioniz_A);
    IonizIntegral Int_B(Int_A);
    EXPECT_TRUE(Int_A == Int_B);

    InterpolationDef InterpolDef;
    IonizInterpolant Interpol_A(Ioniz_A, InterpolDef);
    IonizInterpolant Interpol_B(Interpol_A);
    EXPECT_TRUE(Interpol_A == Interpol_B);
}

// in polymorphism an assignmant and swap operator doesn't make sense

TEST(Ionization, Test_of_dEdx)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Ioniz_dEdx.txt";
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
    double energy;
    double dEdx_stored;
    double dEdx_new;

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> dEdx_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        Ionization Ioniz(particle_def, *medium, ecuts, multiplier);
        IonizIntegral Ioniz_Int(Ioniz);

        dEdx_new = Ioniz_Int.CalculatedEdx(energy);

        ASSERT_NEAR(dEdx_new, dEdx_stored, 1e-10 * dEdx_stored);

        delete medium;
    }
}

TEST(Ionization, Test_of_dNdx)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Ioniz_dNdx.txt";
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
    double energy;
    double dNdx_stored;
    double dNdx_new;

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> dNdx_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        Ionization Ioniz(particle_def, *medium, ecuts, multiplier);
        IonizIntegral Ioniz_Int(Ioniz);

        dNdx_new = Ioniz_Int.CalculatedNdx(energy);

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);

        delete medium;
    }
}

TEST(Ionization, Test_of_dNdx_rnd)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Ioniz_dNdx_rnd.txt";
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
    double energy;
    double rnd;
    double dNdx_rnd_stored;
    double dNdx_rnd_new;

    RandomGenerator::Get().SetSeed(0);

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> rnd >> dNdx_rnd_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        Ionization Ioniz(particle_def, *medium, ecuts, multiplier);
        IonizIntegral Ioniz_Int(Ioniz);

        dNdx_rnd_new = Ioniz_Int.CalculatedNdx(energy, rnd);

        ASSERT_NEAR(dNdx_rnd_new, dNdx_rnd_stored, 1E-10 * dNdx_rnd_stored);

        delete medium;
    }
}

TEST(Ionization, Test_Stochastic_Loss)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Ioniz_e.txt";
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
    double energy;
    double rnd1, rnd2;
    double stochastic_loss_stored;
    double stochastic_loss_new;

    RandomGenerator::Get().SetSeed(0);

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> rnd1 >> rnd2 >>
            stochastic_loss_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        Ionization Ioniz(particle_def, *medium, ecuts, multiplier);
        IonizIntegral Ioniz_Int(Ioniz);

        stochastic_loss_new = Ioniz_Int.CrossSectionIntegral::CalculateStochasticLoss(energy, rnd1, rnd2);

        ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);

        delete medium;
    }
}

TEST(Ionization, Test_of_dEdx_Interpolant)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Ioniz_dEdx_interpol.txt";
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
    double energy;
    double dEdx_stored;
    double dEdx_new;

    InterpolationDef InterpolDef;

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> dEdx_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        Ionization Ioniz(particle_def, *medium, ecuts, multiplier);
        IonizInterpolant Ioniz_Interpol(Ioniz, InterpolDef);

        dEdx_new = Ioniz_Interpol.CalculatedEdx(energy);

        ASSERT_NEAR(dEdx_new, dEdx_stored, 1e-10 * dEdx_stored);

        delete medium;
    }
}

TEST(Ionization, Test_of_dNdx_Interpolant)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Ioniz_dNdx_interpol.txt";
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
    double energy;
    double dNdx_stored;
    double dNdx_new;

    InterpolationDef InterpolDef;

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> dNdx_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        Ionization Ioniz(particle_def, *medium, ecuts, multiplier);
        IonizInterpolant Ioniz_Interpol(Ioniz, InterpolDef);

        dNdx_new = Ioniz_Interpol.CalculatedNdx(energy);

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);

        delete medium;
    }
}

TEST(Ionization, Test_of_dNdxrnd_interpol)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Ioniz_dNdx_rnd_interpol.txt";
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
    double energy;
    double rnd;
    double dNdx_rnd_stored;
    double dNdx_rnd_new;

    InterpolationDef InterpolDef;

    RandomGenerator::Get().SetSeed(0);

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> rnd >> dNdx_rnd_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        Ionization Ioniz(particle_def, *medium, ecuts, multiplier);
        IonizInterpolant Ioniz_Interpol(Ioniz, InterpolDef);

        dNdx_rnd_new = Ioniz_Interpol.CalculatedNdx(energy, rnd);

        ASSERT_NEAR(dNdx_rnd_new, dNdx_rnd_stored, 1E-10 * dNdx_rnd_stored);

        delete medium;
    }
}

TEST(Ionization, Test_of_e_interpol)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Ioniz_e_interpol.txt";
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
    double energy;
    double rnd1, rnd2;
    double stochastic_loss_stored;
    double stochastic_loss_new;

    InterpolationDef InterpolDef;

    RandomGenerator::Get().SetSeed(0);

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> rnd1 >> rnd2 >>
            stochastic_loss_stored;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        Ionization Ioniz(particle_def, *medium, ecuts, multiplier);
        IonizInterpolant Ioniz_Interpol(Ioniz, InterpolDef);

        stochastic_loss_new = Ioniz_Interpol.CrossSectionInterpolant::CalculateStochasticLoss(energy, rnd1, rnd2);

        ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);

        delete medium;
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
