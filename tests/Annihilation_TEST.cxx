
#include "gtest/gtest.h"

#include <fstream>
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/AnnihilationIntegral.h"
#include "PROPOSAL/crossection/AnnihilationInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Annihilation.h"
#include "PROPOSAL/crossection/factories/AnnihilationFactory.h"
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
ParticleDef particle_def = EPlusDef::Get(); //particle
auto medium = std::make_shared<const Water>();
double multiplier   = 1.;

Annihilation* Anni_A = new AnnihilationHeitler(particle_def, medium, multiplier);
Annihilation* Anni_B = new AnnihilationHeitler(particle_def, medium, multiplier);
EXPECT_TRUE(*Anni_A == *Anni_B);

AnnihilationHeitler param_int(particle_def, medium, multiplier);
EXPECT_TRUE(param_int == *Anni_A);

AnnihilationIntegral* Int_A        = new AnnihilationIntegral(param_int);
CrossSectionIntegral* Int_B = new AnnihilationIntegral(param_int);
EXPECT_TRUE(*Int_A == *Int_B);

InterpolationDef InterpolDef;

AnnihilationInterpolant* Interpol_A        = new AnnihilationInterpolant(param_int, InterpolDef);
CrossSectionInterpolant* Interpol_B = new AnnihilationInterpolant(param_int, InterpolDef);
EXPECT_TRUE(*Interpol_A == *Interpol_B);

delete Anni_A;
delete Anni_B;
delete Int_A;
delete Int_B;
delete Interpol_A;
delete Interpol_B;
}

TEST(Comparison, Comparison_not_equal)
{
ParticleDef mu_def  = MuMinusDef::Get();
ParticleDef e_minus_def = EMinusDef::Get();
ParticleDef e_plus_def  = EPlusDef::Get();
auto medium_1 = std::make_shared<const Water>();
auto medium_2 = std::make_shared<const Ice>();
double multiplier_1 = 1.;
double multiplier_2 = 2.;

AnnihilationHeitler Anni_A(mu_def, medium_1, multiplier_1);
AnnihilationHeitler Anni_B(e_minus_def, medium_1, multiplier_1);
AnnihilationHeitler Anni_C(mu_def, medium_2, multiplier_1);
AnnihilationHeitler Anni_D(e_plus_def, medium_1, multiplier_1);
AnnihilationHeitler Anni_E(mu_def, medium_1, multiplier_2);
EXPECT_TRUE(Anni_A != Anni_B);
EXPECT_TRUE(Anni_A != Anni_C);
EXPECT_TRUE(Anni_A != Anni_D);
EXPECT_TRUE(Anni_A != Anni_E);

AnnihilationIntegral Int_A(Anni_A);
AnnihilationIntegral Int_B(Anni_B);
EXPECT_TRUE(Int_A != Int_B);

InterpolationDef InterpolDef;
AnnihilationIntegral AnniIntegral_A(Anni_A);
AnnihilationIntegral AnniIntegral_B(Anni_B);
AnnihilationIntegral AnniIntegral_C(Anni_C);
AnnihilationIntegral AnniIntegral_D(Anni_D);
AnnihilationIntegral AnniIntegral_E(Anni_E);
EXPECT_TRUE(AnniIntegral_A != AnniIntegral_B);
EXPECT_TRUE(AnniIntegral_A != AnniIntegral_C);
EXPECT_TRUE(AnniIntegral_A != AnniIntegral_D);
EXPECT_TRUE(AnniIntegral_A != AnniIntegral_E);

AnnihilationIntegral Integral_A(AnniIntegral_A);
AnnihilationIntegral Integral_B(AnniIntegral_B);
EXPECT_TRUE(Integral_A != Integral_B);
}

TEST(Assignment, Copyconstructor)
{
ParticleDef particle_def = EPlusDef::Get();
auto medium = std::make_shared<const Water>();
double multiplier = 1.;

AnnihilationHeitler Anni_A(particle_def, medium, multiplier);
AnnihilationHeitler Anni_B = Anni_A;
EXPECT_TRUE(Anni_A == Anni_B);

AnnihilationIntegral Int_A(Anni_A);
AnnihilationIntegral Int_B = Int_A;
EXPECT_TRUE(Int_A == Int_B);

InterpolationDef InterpolDef;
AnnihilationInterpolant AnniInterpol_A(Anni_A, InterpolDef);
AnnihilationInterpolant AnniInterpol_B = AnniInterpol_A;
EXPECT_TRUE(AnniInterpol_A == AnniInterpol_B);
}




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
std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);

AnnihilationFactory::Definition annihilation_def;
annihilation_def.multiplier      = multiplier;
annihilation_def.parametrization = AnnihilationFactory::Get().GetEnumFromString(parametrization);

CrossSection* anni = AnnihilationFactory::Get().CreateAnnihilation(particle_def, medium, annihilation_def);
dNdx_new = anni->CalculatedNdx(energy);

ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);

delete anni;
}
}

TEST(Annihilation, Test_of_dNdx_rnd)
{
std::string filename = "bin/TestFiles/Anni_dNdx_rnd.txt";
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
double dNdx_rnd_stored;
double dNdx_rnd_new;

RandomGenerator::Get().SetSeed(0);

while (in.good())
{
in >> particleName >> mediumName >> multiplier >> energy >> parametrization >> rnd >> dNdx_rnd_stored;

ParticleDef particle_def = getParticleDef(particleName);
std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);

AnnihilationFactory::Definition annihilation_def;
annihilation_def.multiplier      = multiplier;
annihilation_def.parametrization = AnnihilationFactory::Get().GetEnumFromString(parametrization);

CrossSection* anni = AnnihilationFactory::Get().CreateAnnihilation(particle_def, medium, annihilation_def);

dNdx_rnd_new = anni->CalculatedNdx(energy, rnd);

ASSERT_NEAR(dNdx_rnd_new, dNdx_rnd_stored, 1E-10 * dNdx_rnd_stored);

delete anni;
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
std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);

AnnihilationFactory::Definition annihilation_def;
annihilation_def.multiplier      = multiplier;
annihilation_def.parametrization = AnnihilationFactory::Get().GetEnumFromString(parametrization);

CrossSection* anni = AnnihilationFactory::Get().CreateAnnihilation(particle_def, medium, annihilation_def);

stochastic_loss_new = anni->CalculateStochasticLoss(energy, rnd1, rnd2);

ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);

delete anni;
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

InterpolationDef InterpolDef;

while (in.good())
{
in >> particleName >> mediumName >> multiplier >> energy >> parametrization >> dNdx_stored;

ParticleDef particle_def = getParticleDef(particleName);
std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);

AnnihilationFactory::Definition annihilation_def;
annihilation_def.multiplier      = multiplier;
annihilation_def.parametrization = AnnihilationFactory::Get().GetEnumFromString(parametrization);

CrossSection* anni = AnnihilationFactory::Get().CreateAnnihilation(particle_def, medium, annihilation_def, InterpolDef);

dNdx_new = anni->CalculatedNdx(energy);

ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);

delete anni;
}
}

TEST(Annihilation, Test_of_dNdxrnd_interpol)
{
std::string filename = "bin/TestFiles/Anni_dNdx_rnd_interpol.txt";
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
double dNdx_rnd_stored;
double dNdx_rnd_new;

InterpolationDef InterpolDef;

RandomGenerator::Get().SetSeed(0);

while (in.good())
{
in >> particleName >> mediumName >> multiplier >>  energy >> parametrization >> rnd >> dNdx_rnd_stored;

ParticleDef particle_def = getParticleDef(particleName);
std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);

AnnihilationFactory::Definition annihilation_def;
annihilation_def.multiplier      = multiplier;
annihilation_def.parametrization = AnnihilationFactory::Get().GetEnumFromString(parametrization);

CrossSection* anni = AnnihilationFactory::Get().CreateAnnihilation(particle_def, medium, annihilation_def, InterpolDef);

dNdx_rnd_new = anni->CalculatedNdx(energy, rnd);

ASSERT_NEAR(dNdx_rnd_new, dNdx_rnd_stored, 1E-10 * dNdx_rnd_stored);

delete anni;
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
std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);

AnnihilationFactory::Definition annihilation_def;
annihilation_def.multiplier      = multiplier;
annihilation_def.parametrization = AnnihilationFactory::Get().GetEnumFromString(parametrization);

CrossSection* anni = AnnihilationFactory::Get().CreateAnnihilation(particle_def, medium, annihilation_def, InterpolDef);

stochastic_loss_new = anni->CalculateStochasticLoss(energy, rnd1, rnd2);

ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);

delete anni;
}
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
