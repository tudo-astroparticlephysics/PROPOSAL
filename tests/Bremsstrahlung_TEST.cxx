
// #include <iostream>
// #include <string>

#include "gtest/gtest.h"

#include <fstream>
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/BremsIntegral.h"
#include "PROPOSAL/crossection/BremsInterpolant.h"
#include "PROPOSAL/crossection/factories/BremsstrahlungFactory.h"
#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"
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

    BremsKelnerKokoulinPetrukhin* Brems_A =
        new BremsKelnerKokoulinPetrukhin(particle_def, medium, ecuts, multiplier, lpm);
    Parametrization* Brems_B = new BremsKelnerKokoulinPetrukhin(particle_def, medium, ecuts, multiplier, lpm);
    EXPECT_TRUE(*Brems_A == *Brems_B);

    BremsKelnerKokoulinPetrukhin param(particle_def, medium, ecuts, multiplier, lpm);
    EXPECT_TRUE(param == *Brems_A);

    BremsIntegral* Int_A        = new BremsIntegral(param);
    CrossSectionIntegral* Int_B = new BremsIntegral(param);
    EXPECT_TRUE(*Int_A == *Int_B);

    InterpolationDef InterpolDef;
    BremsInterpolant* Interpol_A        = new BremsInterpolant(param, InterpolDef);
    CrossSectionInterpolant* Interpol_B = new BremsInterpolant(param, InterpolDef);
    EXPECT_TRUE(*Interpol_A == *Interpol_B);

    delete Brems_A;
    delete Brems_B;
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
    bool lpm_1          = true;
    bool lpm_2          = false;

    BremsKelnerKokoulinPetrukhin Brems_A(mu_def, medium_1, ecuts_1, multiplier_1, lpm_1);
    BremsKelnerKokoulinPetrukhin Brems_B(tau_def, medium_1, ecuts_1, multiplier_1, lpm_1);
    BremsKelnerKokoulinPetrukhin Brems_C(mu_def, medium_2, ecuts_1, multiplier_1, lpm_1);
    BremsKelnerKokoulinPetrukhin Brems_D(mu_def, medium_1, ecuts_2, multiplier_1, lpm_1);
    BremsKelnerKokoulinPetrukhin Brems_E(mu_def, medium_1, ecuts_1, multiplier_2, lpm_1);
    BremsKelnerKokoulinPetrukhin Brems_F(mu_def, medium_1, ecuts_1, multiplier_1, lpm_2);
    EXPECT_TRUE(Brems_A != Brems_B);
    EXPECT_TRUE(Brems_A != Brems_C);
    EXPECT_TRUE(Brems_A != Brems_C);
    EXPECT_TRUE(Brems_A != Brems_E);
    EXPECT_TRUE(Brems_A != Brems_F);

    BremsAndreevBezrukovBugaev param_2(mu_def, medium_1, ecuts_1, multiplier_1, lpm_1);
    BremsPetrukhinShestakov param_3(mu_def, medium_1, ecuts_1, multiplier_1, lpm_1);
    BremsCompleteScreening param_4(mu_def, medium_1, ecuts_1, multiplier_1, lpm_1);
    EXPECT_TRUE(Brems_A != param_2);
    EXPECT_TRUE(Brems_A != param_3);
    EXPECT_TRUE(Brems_A != param_4);
    EXPECT_TRUE(param_2 != param_3);
    EXPECT_TRUE(param_2 != param_4);
    EXPECT_TRUE(param_3 != param_4);

    BremsIntegral Int_A(param_2);
    BremsIntegral Int_B(param_3);
    EXPECT_TRUE(Int_A != Int_B);

    InterpolationDef InterpolDef;
    BremsInterpolant Interpol_A(param_2, InterpolDef);
    BremsInterpolant Interpol_B(param_3, InterpolDef);
    EXPECT_TRUE(Interpol_A != Interpol_B);
}

TEST(Assignment, Copyconstructor)
{
    ParticleDef particle_def = MuMinusDef::Get();
    Water medium;
    EnergyCutSettings ecuts;
    double multiplier = 1.;
    bool lpm          = true;

    BremsKelnerKokoulinPetrukhin Brems_A(particle_def, medium, ecuts, multiplier, lpm);
    BremsKelnerKokoulinPetrukhin Brems_B = Brems_A;
    EXPECT_TRUE(Brems_A == Brems_B);

    BremsIntegral Int_A(Brems_A);
    BremsIntegral Int_B = Int_A;
    EXPECT_TRUE(Int_A == Int_B);

    InterpolationDef InterpolDef;
    BremsInterpolant Interpol_A(Brems_A, InterpolDef);
    BremsInterpolant Interpol_B = Interpol_A;
    EXPECT_TRUE(Interpol_A == Interpol_B);
}

TEST(Assignment, Copyconstructor2)
{
    ParticleDef particle_def = MuMinusDef::Get();
    Water medium;
    EnergyCutSettings ecuts;
    double multiplier = 1.;
    bool lpm          = true;

    BremsKelnerKokoulinPetrukhin Brems_A(particle_def, medium, ecuts, multiplier, lpm);
    BremsKelnerKokoulinPetrukhin Brems_B(Brems_A);
    EXPECT_TRUE(Brems_A == Brems_B);

    BremsIntegral Int_A(Brems_A);
    BremsIntegral Int_B(Int_A);
    EXPECT_TRUE(Int_A == Int_B);

    InterpolationDef InterpolDef;
    BremsInterpolant Interpol_A(Brems_A, InterpolDef);
    BremsInterpolant Interpol_B(Interpol_A);
    EXPECT_TRUE(Interpol_A == Interpol_B);
}

// in polymorphism an assignmant and swap operator doesn't make sense

TEST(Bremsstrahlung, Test_of_dEdx)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Brems_dEdx.txt";
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
    std::string parametrization;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    std::cout.precision(16);

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> dEdx_stored >>
            parametrization;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        BremsstrahlungFactory::Definition brems_def;
        brems_def.multiplier      = multiplier;
        brems_def.lpm_effect      = lpm;
        brems_def.parametrization = BremsstrahlungFactory::Get().GetEnumFromString(parametrization);

        CrossSection* Brems =
            BremsstrahlungFactory::Get().CreateBremsstrahlung(particle_def, *medium, ecuts, brems_def);

        dEdx_new = Brems->CalculatedEdx(energy);

        ASSERT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);

        delete medium;
        delete Brems;
    }
}

TEST(Bremsstrahlung, Test_of_dNdx)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Brems_dNdx.txt";
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
    std::string parametrization;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    std::cout.precision(16);

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> dNdx_stored >>
            parametrization;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        BremsstrahlungFactory::Definition brems_def;
        brems_def.multiplier      = multiplier;
        brems_def.lpm_effect      = lpm;
        brems_def.parametrization = BremsstrahlungFactory::Get().GetEnumFromString(parametrization);

        CrossSection* Brems =
            BremsstrahlungFactory::Get().CreateBremsstrahlung(particle_def, *medium, ecuts, brems_def);

        dNdx_new = Brems->CalculatedNdx(energy);

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);

        delete medium;
        delete Brems;
    }
}

TEST(Bremsstrahlung, Test_of_dNdx_rnd)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Brems_dNdx_rnd.txt";
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
    std::string parametrization;
    double energy;
    double rnd;
    double dNdx_stored;
    double dNdx_new;

    std::cout.precision(16);

    RandomGenerator::Get().SetSeed(0);

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> rnd >> dNdx_stored >>
            parametrization;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        BremsstrahlungFactory::Definition brems_def;
        brems_def.multiplier      = multiplier;
        brems_def.lpm_effect      = lpm;
        brems_def.parametrization = BremsstrahlungFactory::Get().GetEnumFromString(parametrization);

        CrossSection* Brems =
            BremsstrahlungFactory::Get().CreateBremsstrahlung(particle_def, *medium, ecuts, brems_def);

        dNdx_new = Brems->CalculatedNdx(energy, rnd);

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);

        delete medium;
        delete Brems;
    }
}

TEST(Bremsstrahlung, Test_of_e)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Brems_e.txt";
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
    std::string parametrization;
    double energy;
    double rnd1, rnd2;
    double stochastic_loss_stored;
    double stochastic_loss_new;

    std::cout.precision(16);

    RandomGenerator::Get().SetSeed(0);

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> rnd1 >> rnd2 >>
            stochastic_loss_stored >> parametrization;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        BremsstrahlungFactory::Definition brems_def;
        brems_def.multiplier      = multiplier;
        brems_def.lpm_effect      = lpm;
        brems_def.parametrization = BremsstrahlungFactory::Get().GetEnumFromString(parametrization);

        CrossSection* Brems =
            BremsstrahlungFactory::Get().CreateBremsstrahlung(particle_def, *medium, ecuts, brems_def);

        stochastic_loss_new = Brems->CalculateStochasticLoss(energy, rnd1, rnd2);

        ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1e-3 * stochastic_loss_stored);

        delete medium;
        delete Brems;
    }
}

TEST(Bremsstrahlung, Test_of_dEdx_Interpolant)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Brems_dEdx_interpol.txt";
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
    std::string parametrization;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    std::cout.precision(16);
    InterpolationDef InterpolDef;

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> dEdx_stored >>
            parametrization;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        BremsstrahlungFactory::Definition brems_def;
        brems_def.multiplier      = multiplier;
        brems_def.lpm_effect      = lpm;
        brems_def.parametrization = BremsstrahlungFactory::Get().GetEnumFromString(parametrization);

        CrossSection* Brems =
            BremsstrahlungFactory::Get().CreateBremsstrahlung(particle_def, *medium, ecuts, brems_def, InterpolDef);

        dEdx_new = Brems->CalculatedEdx(energy);

        ASSERT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);

        delete medium;
        delete Brems;
    }
}

TEST(Bremsstrahlung, Test_of_dNdx_Interpolant)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Brems_dNdx_interpol.txt";
    in.open(filename.c_str());

    if (!in.good())
    {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    double multiplier;
    bool lpm;
    std::string parametrization;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    std::cout.precision(16);
    InterpolationDef InterpolDef;

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> dNdx_stored >>
            parametrization;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        BremsstrahlungFactory::Definition brems_def;
        brems_def.multiplier      = multiplier;
        brems_def.lpm_effect      = lpm;
        brems_def.parametrization = BremsstrahlungFactory::Get().GetEnumFromString(parametrization);

        CrossSection* Brems =
            BremsstrahlungFactory::Get().CreateBremsstrahlung(particle_def, *medium, ecuts, brems_def, InterpolDef);

        dNdx_new = Brems->CalculatedNdx(energy);

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);

        delete medium;
        delete Brems;
    }
}

TEST(Bremsstrahlung, Test_of_dNdx_rnd_Interpolant)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Brems_dNdx_rnd_interpol.txt";
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
    std::string parametrization;
    double energy;
    double rnd;
    double dNdx_stored;
    double dNdx_new;

    std::cout.precision(16);
    InterpolationDef InterpolDef;

    RandomGenerator::Get().SetSeed(0);

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> rnd >> dNdx_stored >>
            parametrization;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        BremsstrahlungFactory::Definition brems_def;
        brems_def.multiplier      = multiplier;
        brems_def.lpm_effect      = lpm;
        brems_def.parametrization = BremsstrahlungFactory::Get().GetEnumFromString(parametrization);

        CrossSection* Brems =
            BremsstrahlungFactory::Get().CreateBremsstrahlung(particle_def, *medium, ecuts, brems_def, InterpolDef);

        dNdx_new = Brems->CalculatedNdx(energy, rnd);

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);

        delete medium;
        delete Brems;
    }
}

TEST(Bremsstrahlung, Test_of_e_Interpolant)
{
    std::ifstream in;
    std::string filename = testfile_dir + "Brems_e_interpol.txt";
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
    std::string parametrization;
    double energy;
    double rnd1, rnd2;
    double stochastic_loss_stored;
    double stochastic_loss_new;

    std::cout.precision(16);
    InterpolationDef InterpolDef;

    RandomGenerator::Get().SetSeed(0);

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> rnd1 >> rnd2 >>
            stochastic_loss_stored >> parametrization;

        ParticleDef particle_def = getParticleDef(particleName);
        Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings ecuts(ecut, vcut);

        BremsstrahlungFactory::Definition brems_def;
        brems_def.multiplier      = multiplier;
        brems_def.lpm_effect      = lpm;
        brems_def.parametrization = BremsstrahlungFactory::Get().GetEnumFromString(parametrization);

        CrossSection* Brems =
            BremsstrahlungFactory::Get().CreateBremsstrahlung(particle_def, *medium, ecuts, brems_def, InterpolDef);

        stochastic_loss_new = Brems->CalculateStochasticLoss(energy, rnd1, rnd2);

        ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1e-3 * stochastic_loss_stored);

        delete medium;
        delete Brems;
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
