
#include "gtest/gtest.h"

#include <fstream>
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/MupairIntegral.h"
#include "PROPOSAL/crossection/MupairInterpolant.h"
#include "PROPOSAL/crossection/parametrization/MupairProduction.h"
#include "PROPOSAL/crossection/factories/MupairProductionFactory.h"
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
double multiplier   = 1.;

MupairProductionRhoIntegral* MupairInt_A =
        new MupairKelnerKokoulinPetrukhin(particle_def, medium, ecuts, multiplier);
Parametrization* MupairInt_B = new MupairKelnerKokoulinPetrukhin(particle_def, medium, ecuts, multiplier);
EXPECT_TRUE(*MupairInt_A == *MupairInt_B);

MupairKelnerKokoulinPetrukhin param_int(particle_def, medium, ecuts, multiplier);
EXPECT_TRUE(param_int == *MupairInt_A);

MupairIntegral* Int_A        = new MupairIntegral(param_int);
CrossSectionIntegral* Int_B = new MupairIntegral(param_int);
EXPECT_TRUE(*Int_A == *Int_B);

InterpolationDef InterpolDef;
MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin>* MupairInterpol_A =
        new MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin>(particle_def, medium, ecuts, multiplier, InterpolDef);
Parametrization* MupairInterpol_B =
        new MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin>(particle_def, medium, ecuts, multiplier, InterpolDef);
EXPECT_TRUE(*MupairInterpol_A == *MupairInterpol_B);

MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> param_interpol(particle_def, medium, ecuts, multiplier, InterpolDef);
EXPECT_TRUE(param_interpol == *MupairInterpol_A);

MupairInterpolant* Interpol_A        = new MupairInterpolant(param_interpol, InterpolDef);
CrossSectionInterpolant* Interpol_B = new MupairInterpolant(param_interpol, InterpolDef);
EXPECT_TRUE(*Interpol_A == *Interpol_B);

delete MupairInt_A;
delete MupairInt_B;
delete Int_A;
delete Int_B;
delete MupairInterpol_A;
delete MupairInterpol_B;
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

MupairKelnerKokoulinPetrukhin MupairInt_A(mu_def, medium_1, ecuts_1, multiplier_1);
MupairKelnerKokoulinPetrukhin MupairInt_B(tau_def, medium_1, ecuts_1, multiplier_1);
MupairKelnerKokoulinPetrukhin MupairInt_C(mu_def, medium_2, ecuts_1, multiplier_1);
MupairKelnerKokoulinPetrukhin MupairInt_D(mu_def, medium_1, ecuts_2, multiplier_1);
MupairKelnerKokoulinPetrukhin MupairInt_E(mu_def, medium_1, ecuts_1, multiplier_2);
EXPECT_TRUE(MupairInt_A != MupairInt_B);
EXPECT_TRUE(MupairInt_A != MupairInt_C);
EXPECT_TRUE(MupairInt_A != MupairInt_D);
EXPECT_TRUE(MupairInt_A != MupairInt_E);

MupairIntegral Int_A(MupairInt_A);
MupairIntegral Int_B(MupairInt_B);
EXPECT_TRUE(Int_A != Int_B);

InterpolationDef InterpolDef;
MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_A(mu_def, medium_1, ecuts_1, multiplier_1, InterpolDef);
MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_B(tau_def, medium_1, ecuts_1, multiplier_1, InterpolDef);
MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_C(mu_def, medium_2, ecuts_1, multiplier_1, InterpolDef);
MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_D(mu_def, medium_1, ecuts_2, multiplier_1, InterpolDef);
MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_E(mu_def, medium_1, ecuts_1, multiplier_2, InterpolDef);
EXPECT_TRUE(MupairInterpol_A != MupairInterpol_B);
EXPECT_TRUE(MupairInterpol_A != MupairInterpol_C);
EXPECT_TRUE(MupairInterpol_A != MupairInterpol_D);
EXPECT_TRUE(MupairInterpol_A != MupairInterpol_E);

MupairInterpolant Interpol_A(MupairInterpol_A, InterpolDef);
MupairInterpolant Interpol_B(MupairInterpol_B, InterpolDef);
EXPECT_TRUE(Interpol_A != Interpol_B);
}

TEST(Assignment, Copyconstructor)
{
ParticleDef particle_def = MuMinusDef::Get();
Water medium;
EnergyCutSettings ecuts;
double multiplier = 1.;

MupairKelnerKokoulinPetrukhin MupairInt_A(particle_def, medium, ecuts, multiplier);
MupairKelnerKokoulinPetrukhin MupairInt_B = MupairInt_A;
EXPECT_TRUE(MupairInt_A == MupairInt_B);

MupairIntegral Int_A(MupairInt_A);
MupairIntegral Int_B = Int_A;
EXPECT_TRUE(Int_A == Int_B);

InterpolationDef InterpolDef;
MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_A(particle_def, medium, ecuts, multiplier, InterpolDef);
MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_B = MupairInterpol_A;
EXPECT_TRUE(MupairInterpol_A == MupairInterpol_B);

MupairInterpolant Interpol_A(MupairInterpol_A, InterpolDef);
MupairInterpolant Interpol_B = Interpol_A;
EXPECT_TRUE(Interpol_A == Interpol_B);
}

TEST(Assignment, Copyconstructor2)
{
ParticleDef particle_def = MuMinusDef::Get();
Water medium;
EnergyCutSettings ecuts;
double multiplier = 1.;

MupairKelnerKokoulinPetrukhin MupairInt_A(particle_def, medium, ecuts, multiplier);
MupairKelnerKokoulinPetrukhin MupairInt_B(MupairInt_A);
EXPECT_TRUE(MupairInt_A == MupairInt_B);

MupairIntegral Int_A(MupairInt_A);
MupairIntegral Int_B(Int_A);
EXPECT_TRUE(Int_A == Int_B);

InterpolationDef InterpolDef;
MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_A(particle_def, medium, ecuts, multiplier, InterpolDef);
MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_B(MupairInterpol_A);
EXPECT_TRUE(MupairInterpol_A == MupairInterpol_B);

MupairInterpolant Interpol_A(MupairInterpol_A, InterpolDef);
MupairInterpolant Interpol_B(Interpol_A);
EXPECT_TRUE(Interpol_A == Interpol_B);
}

// in polymorphism an assignment and swap operator doesn't make sense

TEST(Mupairproduction, Test_of_dEdx)
{
std::ifstream in;
std::string filename = testfile_dir + "Mupair_dEdx.txt";
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
std::string parametrization;
double dEdx_stored;
double dEdx_new;

while (in.good())
{
in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> parametrization >> dEdx_stored;

ParticleDef particle_def = getParticleDef(particleName);
Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
EnergyCutSettings ecuts(ecut, vcut);

MupairProductionFactory::Definition mupair_def;
mupair_def.multiplier      = multiplier;
mupair_def.parametrization = MupairProductionFactory::Get().GetEnumFromString(parametrization);

CrossSection* mupair = MupairProductionFactory::Get().CreateMupairProduction(particle_def, *medium, ecuts, mupair_def);
dEdx_new = mupair->CalculatedEdx(energy);

ASSERT_NEAR(dEdx_new, dEdx_stored, 1e-10 * dEdx_stored);

delete medium;
delete mupair;
}
}

TEST(Mupairproduction, Test_of_dNdx)
{
std::ifstream in;
std::string filename = testfile_dir + "Mupair_dNdx.txt";
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
std::string parametrization;
double dNdx_stored;
double dNdx_new;

while (in.good())
{
in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> parametrization >> dNdx_stored;

ParticleDef particle_def = getParticleDef(particleName);
Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
EnergyCutSettings ecuts(ecut, vcut);

MupairProductionFactory::Definition mupair_def;
mupair_def.multiplier      = multiplier;
mupair_def.parametrization = MupairProductionFactory::Get().GetEnumFromString(parametrization);

CrossSection* mupair = MupairProductionFactory::Get().CreateMupairProduction(particle_def, *medium, ecuts, mupair_def);
dNdx_new = mupair->CalculatedNdx(energy);

ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);

delete medium;
delete mupair;
}
}

TEST(Mupairproduction, Test_of_dNdx_rnd)
{
std::ifstream in;
std::string filename = testfile_dir + "Mupair_dNdx_rnd.txt";
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
std::string parametrization;
double rnd;
double dNdx_rnd_stored;
double dNdx_rnd_new;

RandomGenerator::Get().SetSeed(0);

while (in.good())
{
in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> parametrization >> rnd >> dNdx_rnd_stored;

ParticleDef particle_def = getParticleDef(particleName);
Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
EnergyCutSettings ecuts(ecut, vcut);

MupairProductionFactory::Definition mupair_def;
mupair_def.multiplier      = multiplier;
mupair_def.parametrization = MupairProductionFactory::Get().GetEnumFromString(parametrization);

CrossSection* mupair = MupairProductionFactory::Get().CreateMupairProduction(particle_def, *medium, ecuts, mupair_def);

dNdx_rnd_new = mupair->CalculatedNdx(energy, rnd);

ASSERT_NEAR(dNdx_rnd_new, dNdx_rnd_stored, 1E-10 * dNdx_rnd_stored);

delete medium;
delete mupair;
}
}

TEST(Mupairproduction, Test_Stochastic_Loss)
{
std::ifstream in;
std::string filename = testfile_dir + "Mupair_e.txt";
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
std::string parametrization;
double rnd1, rnd2;
double stochastic_loss_stored;
double stochastic_loss_new;

std::cout.precision(16);
RandomGenerator::Get().SetSeed(0);

while (in.good())
{
in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> parametrization >> rnd1 >> rnd2 >>
stochastic_loss_stored;


ParticleDef particle_def = getParticleDef(particleName);
Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
EnergyCutSettings ecuts(ecut, vcut);

MupairProductionFactory::Definition mupair_def;
mupair_def.multiplier      = multiplier;
mupair_def.parametrization = MupairProductionFactory::Get().GetEnumFromString(parametrization);

CrossSection* mupair = MupairProductionFactory::Get().CreateMupairProduction(particle_def, *medium, ecuts, mupair_def);

stochastic_loss_new = mupair->CalculateStochasticLoss(energy, rnd1, rnd2);

ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);

delete medium;
delete mupair;
}
}

TEST(Mupairproduction, Test_Calculate_Rho)
{
std::ifstream in;
std::string filename = testfile_dir + "Mupair_rho.txt";
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
double v;
double multiplier;
double energy;
std::string parametrization;
double rnd1, rnd2;
double E1_stored;
double E2_stored;
double E1_new;
double E2_new;

std::vector<Particle*> particles;


std::cout.precision(16);
RandomGenerator::Get().SetSeed(0);

while (in.good())
{
in >> particleName >> mediumName >> ecut >> vcut >> v >> multiplier >> energy >> parametrization >> rnd1 >>
rnd2 >> E1_stored >> E2_stored;


ParticleDef particle_def = getParticleDef(particleName);
Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
EnergyCutSettings ecuts(ecut, vcut);

MupairProductionFactory::Definition mupair_def;
mupair_def.multiplier      = multiplier;
mupair_def.parametrization = MupairProductionFactory::Get().GetEnumFromString(parametrization);

CrossSection* mupair = MupairProductionFactory::Get().CreateMupairProduction(particle_def, *medium, ecuts, mupair_def);

particles = mupair->CalculateProducedParticles(energy, v*energy, rnd1, rnd2);

E1_new = particles[0]->GetEnergy();
E2_new = particles[1]->GetEnergy();
ASSERT_NEAR(E1_new, E1_stored, 1E-6 * E1_stored);
ASSERT_NEAR(E2_new, E2_stored, 1E-6 * E2_stored);

delete medium;
delete mupair;
}
}

TEST(Mupairproduction, Test_of_dEdx_Interpolant)
{
std::ifstream in;
std::string filename = testfile_dir + "Mupair_dEdx_interpol.txt";
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
std::string parametrization;
double dEdx_stored;
double dEdx_new;

InterpolationDef InterpolDef;

while (in.good())
{
in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> parametrization >> dEdx_stored;

ParticleDef particle_def = getParticleDef(particleName);
Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
EnergyCutSettings ecuts(ecut, vcut);

MupairProductionFactory::Definition mupair_def;
mupair_def.multiplier      = multiplier;
mupair_def.parametrization = MupairProductionFactory::Get().GetEnumFromString(parametrization);

CrossSection* mupair = MupairProductionFactory::Get().CreateMupairProduction(particle_def, *medium, ecuts, mupair_def, InterpolDef);

dEdx_new = mupair->CalculatedEdx(energy);

ASSERT_NEAR(dEdx_new, dEdx_stored, 1e-10 * dEdx_stored);

delete medium;
delete mupair;
}
}

TEST(Mupairproduction, Test_of_dNdx_Interpolant)
{
std::ifstream in;
std::string filename = testfile_dir + "Mupair_dNdx_interpol.txt";
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
std::string parametrization;
double dNdx_stored;
double dNdx_new;

InterpolationDef InterpolDef;

while (in.good())
{
in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> parametrization >> dNdx_stored;

ParticleDef particle_def = getParticleDef(particleName);
Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
EnergyCutSettings ecuts(ecut, vcut);

MupairProductionFactory::Definition mupair_def;
mupair_def.multiplier      = multiplier;
mupair_def.parametrization = MupairProductionFactory::Get().GetEnumFromString(parametrization);

CrossSection* mupair = MupairProductionFactory::Get().CreateMupairProduction(particle_def, *medium, ecuts, mupair_def, InterpolDef);

dNdx_new = mupair->CalculatedNdx(energy);

ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);

delete medium;
delete mupair;
}
}

TEST(Mupairproduction, Test_of_dNdxrnd_interpol)
{
std::ifstream in;
std::string filename = testfile_dir + "Mupair_dNdx_rnd_interpol.txt";
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
std::string parametrization;
double rnd;
double dNdx_rnd_stored;
double dNdx_rnd_new;

InterpolationDef InterpolDef;

RandomGenerator::Get().SetSeed(0);

while (in.good())
{
in >> particleName >> mediumName >> ecut >> vcut >> multiplier >>  energy >> parametrization >> rnd >> dNdx_rnd_stored;

ParticleDef particle_def = getParticleDef(particleName);
Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
EnergyCutSettings ecuts(ecut, vcut);

MupairProductionFactory::Definition mupair_def;
mupair_def.multiplier      = multiplier;
mupair_def.parametrization = MupairProductionFactory::Get().GetEnumFromString(parametrization);

CrossSection* mupair = MupairProductionFactory::Get().CreateMupairProduction(particle_def, *medium, ecuts, mupair_def, InterpolDef);

dNdx_rnd_new = mupair->CalculatedNdx(energy, rnd);

ASSERT_NEAR(dNdx_rnd_new, dNdx_rnd_stored, 1E-10 * dNdx_rnd_stored);

delete medium;
delete mupair;
}
}

TEST(Mupairproduction, Test_of_e_interpol)
{
std::ifstream in;
std::string filename = testfile_dir + "Mupair_e_interpol.txt";
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
std::string parametrization;
double rnd1, rnd2;
double stochastic_loss_stored;
double stochastic_loss_new;

InterpolationDef InterpolDef;

RandomGenerator::Get().SetSeed(0);

while (in.good())
{
in >> particleName >> mediumName >> ecut >> vcut >> multiplier >>  energy >> parametrization >> rnd1 >> rnd2 >>
stochastic_loss_stored;

ParticleDef particle_def = getParticleDef(particleName);
Medium* medium           = MediumFactory::Get().CreateMedium(mediumName);
EnergyCutSettings ecuts(ecut, vcut);

MupairProductionFactory::Definition mupair_def;
mupair_def.multiplier      = multiplier;
mupair_def.parametrization = MupairProductionFactory::Get().GetEnumFromString(parametrization);

CrossSection* mupair = MupairProductionFactory::Get().CreateMupairProduction(particle_def, *medium, ecuts, mupair_def, InterpolDef);

stochastic_loss_new = mupair->CalculateStochasticLoss(energy, rnd1, rnd2);

ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);

delete medium;
delete mupair;
}
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
