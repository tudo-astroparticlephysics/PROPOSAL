
#include "gtest/gtest.h"

#include <fstream>
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/WeakIntegral.h"
#include "PROPOSAL/crossection/WeakInterpolant.h"
#include "PROPOSAL/crossection/parametrization/WeakInteraction.h"
#include "PROPOSAL/crossection/factories/WeakInteractionFactory.h"
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
    } else if (name == "EMinus")
    {
        return EMinusDef::Get();
    } else if (name == "MuPlus")
    {
        return MuPlusDef::Get();
    } else if (name == "TauPlus")
    {
        return TauPlusDef::Get();
    } else if (name == "EPlus")
    {
        return EPlusDef::Get();
    }
    else{
        return MuMinusDef::Get();
    }
}

TEST(Comparison, Comparison_equal_particle)
{
ParticleDef particle_def = MuMinusDef::Get(); //particle
Water medium;
EnergyCutSettings ecuts;
double multiplier   = 1.;

WeakInteraction* Weak_A = new WeakCooperSarkarMertsch(particle_def, medium, multiplier);
Parametrization* Weak_B = new WeakCooperSarkarMertsch(particle_def, medium, multiplier);
EXPECT_TRUE(*Weak_B == *Weak_B);

WeakCooperSarkarMertsch param_int(particle_def, medium, multiplier);
EXPECT_TRUE(param_int == *Weak_A);

WeakIntegral* Int_A        = new WeakIntegral(param_int);
CrossSectionIntegral* Int_B = new WeakIntegral(param_int);
EXPECT_TRUE(*Int_A == *Int_B);

InterpolationDef InterpolDef;

WeakInterpolant* Interpol_A        = new WeakInterpolant(param_int, InterpolDef);
CrossSectionInterpolant* Interpol_B = new WeakInterpolant(param_int, InterpolDef);
EXPECT_TRUE(*Interpol_A == *Interpol_B);

delete Weak_A;
delete Weak_B;
delete Int_A;
delete Int_B;
delete Interpol_A;
delete Interpol_B;
}

TEST(Comparison, Comparison_equal_antiparticle)
{
    ParticleDef particle_def = MuPlusDef::Get(); //antiparticle
    Water medium;
    double multiplier   = 1.;

    WeakInteraction* Weak_A = new WeakCooperSarkarMertsch(particle_def, medium, multiplier);
    Parametrization* Weak_B = new WeakCooperSarkarMertsch(particle_def, medium, multiplier);
    EXPECT_TRUE(*Weak_B == *Weak_B);

    WeakCooperSarkarMertsch param_int(particle_def, medium, multiplier);
    EXPECT_TRUE(param_int == *Weak_A);

    WeakIntegral* Int_A        = new WeakIntegral(param_int);
    CrossSectionIntegral* Int_B = new WeakIntegral(param_int);
    EXPECT_TRUE(*Int_A == *Int_B);

    InterpolationDef InterpolDef;

    WeakInterpolant* Interpol_A        = new WeakInterpolant(param_int, InterpolDef);
    CrossSectionInterpolant* Interpol_B = new WeakInterpolant(param_int, InterpolDef);
    EXPECT_TRUE(*Interpol_A == *Interpol_B);

    delete Weak_A;
    delete Weak_B;
    delete Int_A;
    delete Int_B;
    delete Interpol_A;
    delete Interpol_B;
}

TEST(Comparison, Comparison_not_equal)
{
ParticleDef mu_def  = MuMinusDef::Get();
ParticleDef tau_def = TauMinusDef::Get();
ParticleDef mu_plus_def  = MuPlusDef::Get();
Water medium_1;
Ice medium_2;
double multiplier_1 = 1.;
double multiplier_2 = 2.;

WeakCooperSarkarMertsch Weak_A(mu_def, medium_1, multiplier_1);
WeakCooperSarkarMertsch Weak_B(tau_def, medium_1, multiplier_1);
WeakCooperSarkarMertsch Weak_C(mu_def, medium_2, multiplier_1);
WeakCooperSarkarMertsch Weak_D(mu_plus_def, medium_1, multiplier_1);
WeakCooperSarkarMertsch Weak_E(mu_def, medium_1, multiplier_2);
EXPECT_TRUE(Weak_A != Weak_B);
EXPECT_TRUE(Weak_A != Weak_C);
EXPECT_TRUE(Weak_A != Weak_D);
EXPECT_TRUE(Weak_A != Weak_E);

WeakIntegral Int_A(Weak_A);
WeakIntegral Int_B(Weak_B);
EXPECT_TRUE(Int_A != Int_B);

InterpolationDef InterpolDef;
WeakIntegral WeakIntegral_A(Weak_A);
WeakIntegral WeakIntegral_B(Weak_B);
WeakIntegral WeakIntegral_C(Weak_C);
WeakIntegral WeakIntegral_D(Weak_D);
WeakIntegral WeakIntegral_E(Weak_E);
EXPECT_TRUE(WeakIntegral_A != WeakIntegral_B);
EXPECT_TRUE(WeakIntegral_A != WeakIntegral_C);
EXPECT_TRUE(WeakIntegral_A != WeakIntegral_D);
EXPECT_TRUE(WeakIntegral_A != WeakIntegral_E);

WeakIntegral Integral_A(WeakIntegral_A);
WeakIntegral Integral_B(WeakIntegral_B);
EXPECT_TRUE(Integral_A != Integral_B);
}

TEST(Assignment, Copyconstructor)
{
ParticleDef particle_def = MuMinusDef::Get();
Water medium;
double multiplier = 1.;

WeakCooperSarkarMertsch Weak_A(particle_def, medium, multiplier);
WeakCooperSarkarMertsch Weak_B = Weak_A;
EXPECT_TRUE(Weak_A == Weak_B);

WeakIntegral Int_A(Weak_A);
WeakIntegral Int_B = Int_A;
EXPECT_TRUE(Int_A == Int_B);

InterpolationDef InterpolDef;
WeakInterpolant WeakInterpol_A(Weak_A, InterpolDef);
WeakInterpolant WeakInterpol_B = WeakInterpol_A;
EXPECT_TRUE(WeakInterpol_A == WeakInterpol_B);
}




TEST(WeakInteraction, Test_of_dNdx)
{
std::ifstream in;
std::string filename = "bin/TestFiles/Weak_dNdx.txt";
in.open(filename.c_str());

if (!in.good())
{
std::cerr << "File \"" << filename << "\" not found" << std::endl;
}

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
Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);

WeakInteractionFactory::Definition weak_def;
weak_def.multiplier      = multiplier;
weak_def.parametrization = WeakInteractionFactory::Get().GetEnumFromString(parametrization);

CrossSection* weak = WeakInteractionFactory::Get().CreateWeakInteraction(particle_def, *medium, weak_def);
dNdx_new = weak->CalculatedNdx(energy);

ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);

delete medium;
delete weak;
}
}

TEST(WeakInteraction, Test_of_dNdx_rnd)
{
std::ifstream in;
std::string filename = "bin/TestFiles/Weak_dNdx_rnd.txt";
in.open(filename.c_str());

if (!in.good())
{
std::cerr << "File \"" << filename << "\" not found" << std::endl;
}

char firstLine[256];
in.getline(firstLine, 256);

std::string particleName;
std::string mediumName;
double multiplier;
double energy;
std::string parametrization;
double rnd;
double dNdx_rnd_stored;
double dNdx_rnd_new;

RandomGenerator::Get().SetSeed(0);

while (in.good())
{
in >> particleName >> mediumName >> multiplier >> energy >> parametrization >> rnd >> dNdx_rnd_stored;

ParticleDef particle_def = getParticleDef(particleName);
Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);

WeakInteractionFactory::Definition weak_def;
weak_def.multiplier      = multiplier;
weak_def.parametrization = WeakInteractionFactory::Get().GetEnumFromString(parametrization);

CrossSection* weak = WeakInteractionFactory::Get().CreateWeakInteraction(particle_def, *medium, weak_def);

dNdx_rnd_new = weak->CalculatedNdx(energy, rnd);

ASSERT_NEAR(dNdx_rnd_new, dNdx_rnd_stored, 1E-10 * dNdx_rnd_stored);

delete medium;
delete weak;
}
}

TEST(WeakInteraction, Test_Stochastic_Loss)
{
std::ifstream in;
std::string filename = "bin/TestFiles/Weak_e.txt";
in.open(filename.c_str());

if (!in.good())
{
std::cerr << "File \"" << filename << "\" not found" << std::endl;
}

char firstLine[256];
in.getline(firstLine, 256);

std::string particleName;
std::string mediumName;
double multiplier;
double energy;
std::string parametrization;
double rnd1, rnd2;
double stochastic_loss_stored;
double stochastic_loss_new;

std::cout.precision(16);
RandomGenerator::Get().SetSeed(0);

while (in.good())
{
in >> particleName >> mediumName >> multiplier >> energy >> parametrization >> rnd1 >> rnd2 >>
stochastic_loss_stored;


ParticleDef particle_def = getParticleDef(particleName);
Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);

WeakInteractionFactory::Definition weak_def;
weak_def.multiplier      = multiplier;
weak_def.parametrization = WeakInteractionFactory::Get().GetEnumFromString(parametrization);

CrossSection* weak = WeakInteractionFactory::Get().CreateWeakInteraction(particle_def, *medium, weak_def);

stochastic_loss_new = weak->CalculateStochasticLoss(energy, rnd1, rnd2);

ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);

delete medium;
delete weak;
}
}


TEST(WeakInteraction, Test_of_dNdx_Interpolant)
{
std::ifstream in;
std::string filename = "bin/TestFiles/Weak_dNdx_interpol.txt";
in.open(filename.c_str());

if (!in.good())
{
std::cerr << "File \"" << filename << "\" not found" << std::endl;
}

char firstLine[256];
in.getline(firstLine, 256);

std::string particleName;
std::string mediumName;
double multiplier;
double energy;
std::string parametrization;
double dNdx_stored;
double dNdx_new;

InterpolationDef InterpolDef;

while (in.good())
{
in >> particleName >> mediumName >> multiplier >> energy >> parametrization >> dNdx_stored;

ParticleDef particle_def = getParticleDef(particleName);
Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);

WeakInteractionFactory::Definition weak_def;
weak_def.multiplier      = multiplier;
weak_def.parametrization = WeakInteractionFactory::Get().GetEnumFromString(parametrization);

CrossSection* weak = WeakInteractionFactory::Get().CreateWeakInteraction(particle_def, *medium, weak_def, InterpolDef);

dNdx_new = weak->CalculatedNdx(energy);

ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);

delete medium;
delete weak;
}
}

TEST(WeakInteraction, Test_of_dNdxrnd_interpol)
{
std::ifstream in;
std::string filename = "bin/TestFiles/Weak_dNdx_rnd_interpol.txt";
in.open(filename.c_str());

if (!in.good())
{
std::cerr << "File \"" << filename << "\" not found" << std::endl;
}

char firstLine[256];
in.getline(firstLine, 256);

std::string particleName;
std::string mediumName;
double multiplier;
double energy;
std::string parametrization;
double rnd;
double dNdx_rnd_stored;
double dNdx_rnd_new;

InterpolationDef InterpolDef;

RandomGenerator::Get().SetSeed(0);

while (in.good())
{
in >> particleName >> mediumName >> multiplier >>  energy >> parametrization >> rnd >> dNdx_rnd_stored;

ParticleDef particle_def = getParticleDef(particleName);
Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);

WeakInteractionFactory::Definition weak_def;
weak_def.multiplier      = multiplier;
weak_def.parametrization = WeakInteractionFactory::Get().GetEnumFromString(parametrization);

CrossSection* weak = WeakInteractionFactory::Get().CreateWeakInteraction(particle_def, *medium, weak_def, InterpolDef);

dNdx_rnd_new = weak->CalculatedNdx(energy, rnd);

ASSERT_NEAR(dNdx_rnd_new, dNdx_rnd_stored, 1E-10 * dNdx_rnd_stored);

delete medium;
delete weak;
}
}

TEST(WeakInteraction, Test_of_e_interpol)
{
std::ifstream in;
std::string filename = "bin/TestFiles/Weak_e_interpol.txt";
in.open(filename.c_str());

if (!in.good())
{
std::cerr << "File \"" << filename << "\" not found" << std::endl;
}

char firstLine[256];
in.getline(firstLine, 256);

std::string particleName;
std::string mediumName;
double multiplier;
double energy;
std::string parametrization;
double rnd1, rnd2;
double stochastic_loss_stored;
double stochastic_loss_new;

InterpolationDef InterpolDef;

RandomGenerator::Get().SetSeed(0);

while (in.good())
{
in >> particleName >> mediumName >> multiplier >>  energy >> parametrization >> rnd1 >> rnd2 >>
stochastic_loss_stored;

ParticleDef particle_def = getParticleDef(particleName);
Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);

WeakInteractionFactory::Definition weak_def;
weak_def.multiplier      = multiplier;
weak_def.parametrization = WeakInteractionFactory::Get().GetEnumFromString(parametrization);

CrossSection* weak = WeakInteractionFactory::Get().CreateWeakInteraction(particle_def, *medium, weak_def, InterpolDef);

stochastic_loss_new = weak->CalculateStochasticLoss(energy, rnd1, rnd2);

ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);

delete medium;
delete weak;
}
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
