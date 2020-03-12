
// #include <iostream>
// #include <string>

#include "gtest/gtest.h"

#include <fstream>
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/ComptonIntegral.h"
#include "PROPOSAL/crossection/ComptonInterpolant.h"
#include "PROPOSAL/crossection/factories/ComptonFactory.h"
#include "PROPOSAL/crossection/parametrization/Compton.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

const std::string testfile_dir = "bin/TestFiles/";

TEST(Comparison, Comparison_equal)
{
ParticleDef particle_def = GammaDef::Get();
auto medium = std::make_shared<const Water>();
EnergyCutSettings ecuts;
double multiplier = 1.;

ComptonKleinNishina* Compton_A =
        new ComptonKleinNishina(particle_def, medium, ecuts, multiplier);
Parametrization* Compton_B = new ComptonKleinNishina(particle_def, medium, ecuts, multiplier);
EXPECT_TRUE(*Compton_A == *Compton_B);

ComptonKleinNishina param(particle_def, medium, ecuts, multiplier);
EXPECT_TRUE(param == *Compton_A);

ComptonIntegral* Int_A        = new ComptonIntegral(param);
CrossSectionIntegral* Int_B = new ComptonIntegral(param);
EXPECT_TRUE(*Int_A == *Int_B);

InterpolationDef InterpolDef;
ComptonInterpolant* Interpol_A        = new ComptonInterpolant(param, InterpolDef);
CrossSectionInterpolant* Interpol_B = new ComptonInterpolant(param, InterpolDef);
EXPECT_TRUE(*Interpol_A == *Interpol_B);

delete Compton_A;
delete Compton_B;
delete Int_A;
delete Int_B;
delete Interpol_A;
delete Interpol_B;
}

TEST(Comparison, Comparison_not_equal)
{
ParticleDef gamma_def  = GammaDef::Get();
auto medium_1 = std::make_shared<const Water>();
auto medium_2 = std::make_shared<const Ice>();
EnergyCutSettings ecuts_1(500, -1);
EnergyCutSettings ecuts_2(-1, 0.05);
double multiplier_1 = 1.;
double multiplier_2 = 2.;

ComptonKleinNishina Compton_A(gamma_def, medium_1, ecuts_1, multiplier_1);
ComptonKleinNishina Compton_B(gamma_def, medium_2, ecuts_1, multiplier_1);
ComptonKleinNishina Compton_C(gamma_def, medium_1, ecuts_2, multiplier_1);
ComptonKleinNishina Compton_D(gamma_def, medium_1, ecuts_1, multiplier_2);
EXPECT_TRUE(Compton_A != Compton_B);
EXPECT_TRUE(Compton_A != Compton_C);
EXPECT_TRUE(Compton_A != Compton_C);
}

TEST(Assignment, Copyconstructor)
{
ParticleDef particle_def = GammaDef::Get();
auto medium = std::make_shared<const Water>();
EnergyCutSettings ecuts;
double multiplier = 1.;

ComptonKleinNishina Compton_A(particle_def, medium, ecuts, multiplier);
ComptonKleinNishina Compton_B = Compton_A;
EXPECT_TRUE(Compton_A == Compton_B);

ComptonIntegral Int_A(Compton_A);
ComptonIntegral Int_B = Int_A;
EXPECT_TRUE(Int_A == Int_B);

InterpolationDef InterpolDef;
ComptonInterpolant Interpol_A(Compton_A, InterpolDef);
ComptonInterpolant Interpol_B = Interpol_A;
EXPECT_TRUE(Interpol_A == Interpol_B);
}

TEST(Assignment, Copyconstructor2)
{
ParticleDef particle_def = GammaDef::Get();
auto medium = std::make_shared<const Water>();
EnergyCutSettings ecuts;
double multiplier = 1.;

ComptonKleinNishina Compton_A(particle_def, medium, ecuts, multiplier);
ComptonKleinNishina Compton_B(Compton_A);
EXPECT_TRUE(Compton_A == Compton_B);

ComptonIntegral Int_A(Compton_A);
ComptonIntegral Int_B(Int_A);
EXPECT_TRUE(Int_A == Int_B);

InterpolationDef InterpolDef;
ComptonInterpolant Interpol_A(Compton_A, InterpolDef);
ComptonInterpolant Interpol_B(Interpol_A);
EXPECT_TRUE(Interpol_A == Interpol_B);
}

// in polymorphism an assignment and swap operator doesn't make sense

TEST(Compton, Test_of_dEdx)
{
std::string filename = testfile_dir + "Compton_dEdx.txt";
std::ifstream in{filename};
EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

char firstLine[256];
in.getline(firstLine, 256);

std::string mediumName;
double ecut;
double vcut;
double multiplier;
std::string parametrization;
double energy;
double dEdx_stored;
double dEdx_new;

std::cout.precision(16);

while (in.good())
{
in >> mediumName >> ecut >> vcut >> multiplier >> energy >> dEdx_stored >>
parametrization;

ParticleDef particle_def = GammaDef::Get();
std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);
EnergyCutSettings ecuts(ecut, vcut);

ComptonFactory::Definition compton_def;
compton_def.multiplier      = multiplier;
compton_def.parametrization = ComptonFactory::Get().GetEnumFromString(parametrization);

CrossSection* Compton =
        ComptonFactory::Get().CreateCompton(particle_def, medium, ecuts, compton_def);

dEdx_new = Compton->CalculatedEdx(energy);

ASSERT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);

delete Compton;
}
}

TEST(Compton, Test_of_dNdx)
{
std::string filename = testfile_dir + "Compton_dNdx.txt";
std::ifstream in{filename};
EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

char firstLine[256];
in.getline(firstLine, 256);

std::string mediumName;
double ecut;
double vcut;
double multiplier;
std::string parametrization;
double energy;
double dNdx_stored;
double dNdx_new;

std::cout.precision(16);

while (in.good())
{
in >> mediumName >> ecut >> vcut >> multiplier >> energy >> dNdx_stored >>
parametrization;

ParticleDef particle_def = GammaDef::Get();
std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);
EnergyCutSettings ecuts(ecut, vcut);

ComptonFactory::Definition compton_def;
compton_def.multiplier      = multiplier;
compton_def.parametrization = ComptonFactory::Get().GetEnumFromString(parametrization);

CrossSection* Compton =
        ComptonFactory::Get().CreateCompton(particle_def, medium, ecuts, compton_def);

dNdx_new = Compton->CalculatedNdx(energy);

ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);

delete Compton;
}
}

TEST(Compton, Test_of_dNdx_rnd)
{
std::string filename = testfile_dir + "Compton_dNdx_rnd.txt";
std::ifstream in{filename};
EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

char firstLine[256];
in.getline(firstLine, 256);

std::string mediumName;
double ecut;
double vcut;
double multiplier;
std::string parametrization;
double energy;
double rnd;
double dNdx_stored;
double dNdx_new;

std::cout.precision(16);

RandomGenerator::Get().SetSeed(0);

while (in.good())
{
in >> mediumName >> ecut >> vcut >> multiplier >> energy >> rnd >> dNdx_stored >>
parametrization;

ParticleDef particle_def = GammaDef::Get();
std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);
EnergyCutSettings ecuts(ecut, vcut);

ComptonFactory::Definition compton_def;
compton_def.multiplier      = multiplier;
compton_def.parametrization = ComptonFactory::Get().GetEnumFromString(parametrization);

CrossSection* Compton =
        ComptonFactory::Get().CreateCompton(particle_def, medium, ecuts, compton_def);

dNdx_new = Compton->CalculatedNdx(energy, rnd);

ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);

delete Compton;
}
}

TEST(Compton, Test_of_e)
{
std::string filename = testfile_dir + "Compton_e.txt";
std::ifstream in{filename};
EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

char firstLine[256];
in.getline(firstLine, 256);

std::string mediumName;
double ecut;
double vcut;
double multiplier;
std::string parametrization;
double energy;
double rnd1, rnd2;
double stochastic_loss_stored;
double stochastic_loss_new;

std::cout.precision(16);

RandomGenerator::Get().SetSeed(0);

while (in.good())
{
in >> mediumName >> ecut >> vcut >> multiplier >> energy >> rnd1 >> rnd2 >>
stochastic_loss_stored >> parametrization;

ParticleDef particle_def = GammaDef::Get();
std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);
EnergyCutSettings ecuts(ecut, vcut);

ComptonFactory::Definition compton_def;
compton_def.multiplier      = multiplier;
compton_def.parametrization = ComptonFactory::Get().GetEnumFromString(parametrization);

CrossSection* Compton =
        ComptonFactory::Get().CreateCompton(particle_def, medium, ecuts, compton_def);

stochastic_loss_new = Compton->CalculateStochasticLoss(energy, rnd1, rnd2);

ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1e-3 * stochastic_loss_stored);

delete Compton;
}
}

TEST(Compton, Test_of_dEdx_Interpolant)
{
std::string filename = testfile_dir + "Compton_dEdx_interpol.txt";
std::ifstream in{filename};
EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

char firstLine[256];
in.getline(firstLine, 256);

std::string mediumName;
double ecut;
double vcut;
double multiplier;
std::string parametrization;
double energy;
double dEdx_stored;
double dEdx_new;

std::cout.precision(16);
InterpolationDef InterpolDef;

while (in.good())
{
in >> mediumName >> ecut >> vcut >> multiplier >> energy >> dEdx_stored >>
parametrization;

ParticleDef particle_def = GammaDef::Get();
std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);
EnergyCutSettings ecuts(ecut, vcut);

ComptonFactory::Definition compton_def;
compton_def.multiplier      = multiplier;
compton_def.parametrization = ComptonFactory::Get().GetEnumFromString(parametrization);

CrossSection* Compton =
        ComptonFactory::Get().CreateCompton(particle_def, medium, ecuts, compton_def, InterpolDef);

dEdx_new = Compton->CalculatedEdx(energy);

ASSERT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);

delete Compton;
}
}

TEST(Compton, Test_of_dNdx_Interpolant)
{
std::string filename = testfile_dir + "Compton_dNdx_interpol.txt";
std::ifstream in{filename};
EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

std::string mediumName;
double ecut;
double vcut;
double multiplier;
std::string parametrization;
double energy;
double dNdx_stored;
double dNdx_new;

std::cout.precision(16);
InterpolationDef InterpolDef;

while (in.good())
{
in >> mediumName >> ecut >> vcut >> multiplier >> energy >> dNdx_stored >> parametrization;

ParticleDef particle_def = GammaDef::Get();
std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);
EnergyCutSettings ecuts(ecut, vcut);

ComptonFactory::Definition compton_def;
compton_def.multiplier      = multiplier;
compton_def.parametrization = ComptonFactory::Get().GetEnumFromString(parametrization);

CrossSection* Compton =
        ComptonFactory::Get().CreateCompton(particle_def, medium, ecuts, compton_def, InterpolDef);

dNdx_new = Compton->CalculatedNdx(energy);

ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);

delete Compton;
}
}

TEST(Compton, Test_of_dNdx_rnd_Interpolant)
{
std::string filename = testfile_dir + "Compton_dNdx_rnd_interpol.txt";
std::ifstream in{filename};
EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

char firstLine[256];
in.getline(firstLine, 256);

std::string mediumName;
double ecut;
double vcut;
double multiplier;
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
in >> mediumName >> ecut >> vcut >> multiplier >> energy >> rnd >> dNdx_stored >> parametrization;

ParticleDef particle_def = GammaDef::Get();
std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);
EnergyCutSettings ecuts(ecut, vcut);

ComptonFactory::Definition compton_def;
compton_def.multiplier      = multiplier;
compton_def.parametrization = ComptonFactory::Get().GetEnumFromString(parametrization);

CrossSection* Compton =
        ComptonFactory::Get().CreateCompton(particle_def, medium, ecuts, compton_def, InterpolDef);

dNdx_new = Compton->CalculatedNdx(energy, rnd);

ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);

delete Compton;
}
}

TEST(Compton, Test_of_e_Interpolant)
{
std::string filename = testfile_dir + "Compton_e_interpol.txt";
std::ifstream in{filename};
EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

char firstLine[256];
in.getline(firstLine, 256);

std::string mediumName;
double ecut;
double vcut;
double multiplier;
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
in >> mediumName >> ecut >> vcut >> multiplier >> energy >> rnd1 >> rnd2 >> stochastic_loss_stored >> parametrization;

ParticleDef particle_def = GammaDef::Get();
std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);
EnergyCutSettings ecuts(ecut, vcut);

ComptonFactory::Definition compton_def;
compton_def.multiplier      = multiplier;
compton_def.parametrization = ComptonFactory::Get().GetEnumFromString(parametrization);

CrossSection* Compton =
        ComptonFactory::Get().CreateCompton(particle_def, medium, ecuts, compton_def, InterpolDef);

stochastic_loss_new = Compton->CalculateStochasticLoss(energy, rnd1, rnd2);

ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1e-3 * stochastic_loss_stored);

delete Compton;
}
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
