
#include "gtest/gtest.h"

#include <fstream>
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/PhotoPairIntegral.h"
#include "PROPOSAL/crossection/PhotoPairInterpolant.h"
#include "PROPOSAL/crossection/parametrization/PhotoPairProduction.h"
#include "PROPOSAL/crossection/factories/PhotoPairFactory.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "Gamma")
    {
        return GammaDef::Get();
    }
    else{
        EXPECT_TRUE(false);
        return MuMinusDef::Get();
    }
}

TEST(Comparison, Comparison_equal)
{
ParticleDef particle_def = GammaDef::Get();
auto medium = std::make_shared<const Water>();
double multiplier   = 1.;

PhotoPairProduction* PhotoPair_A = new PhotoPairTsai(particle_def, medium, multiplier);
Parametrization* PhotoPair_B = new PhotoPairTsai(particle_def, medium, multiplier);
EXPECT_TRUE(*PhotoPair_A == *PhotoPair_B);

PhotoPairTsai param_int(particle_def, medium, multiplier);
EXPECT_TRUE(param_int == *PhotoPair_A);

PhotoAngleNoDeflection const PhotoAngle(particle_def, medium);

PhotoPairIntegral* Int_A        = new PhotoPairIntegral(param_int, PhotoAngle);
CrossSectionIntegral* Int_B = new PhotoPairIntegral(param_int, PhotoAngle);
EXPECT_TRUE(*Int_A == *Int_B);

InterpolationDef InterpolDef;

PhotoPairInterpolant* Interpol_A        = new PhotoPairInterpolant(param_int, PhotoAngle, InterpolDef);
CrossSectionInterpolant* Interpol_B = new PhotoPairInterpolant(param_int, PhotoAngle, InterpolDef);
EXPECT_TRUE(*Interpol_A == *Interpol_B);

delete PhotoPair_A;
delete PhotoPair_B;
delete Int_A;
delete Int_B;
delete Interpol_A;
delete Interpol_B;
}

TEST(Comparison, Comparison_equal_PhotoAngleDistribution)
{
ParticleDef particle_def = GammaDef::Get(); //particle
auto medium = std::make_shared<const Water>();

PhotoAngleDistribution* PhotoAngle_A = new PhotoAngleNoDeflection(particle_def, medium);
PhotoAngleNoDeflection* PhotoAngle_B = new PhotoAngleNoDeflection(particle_def, medium);

EXPECT_TRUE(*PhotoAngle_A == *PhotoAngle_B);

PhotoAngleNoDeflection angle(particle_def, medium);
EXPECT_TRUE(angle == *PhotoAngle_A);

delete PhotoAngle_A;
delete PhotoAngle_B;
}

TEST(Comparison, Comparison_not_equal)
{
ParticleDef photon  = GammaDef::Get();
ParticleDef mu_def  = MuMinusDef::Get(); //makes no physical sense but programmatically possible
auto medium_1 = std::make_shared<const Water>();
auto medium_2 = std::make_shared<const Ice>();
double multiplier_1 = 1.;
double multiplier_2 = 2.;
PhotoAngleNoDeflection const PhotoAngle_1(photon, medium_1);
PhotoAngleEGS const PhotoAngle_2(photon, medium_1);

PhotoPairTsai PhotoPair_A(photon, medium_1, multiplier_1);
PhotoPairTsai PhotoPair_B(photon, medium_2, multiplier_1);
PhotoPairTsai PhotoPair_C(photon, medium_1, multiplier_2);
PhotoPairTsai PhotoPair_D(mu_def, medium_1, multiplier_1);
EXPECT_TRUE(PhotoPair_A != PhotoPair_B);
EXPECT_TRUE(PhotoPair_A != PhotoPair_C);
EXPECT_TRUE(PhotoPair_A != PhotoPair_D);

PhotoPairIntegral Int_A(PhotoPair_A, PhotoAngle_1);
PhotoPairIntegral Int_B(PhotoPair_A, PhotoAngle_2);
PhotoPairIntegral Int_C(PhotoPair_B, PhotoAngle_1);
EXPECT_TRUE(Int_A != Int_B);
EXPECT_TRUE(Int_A != Int_C);

PhotoPairIntegral Integral_A(Int_A);
PhotoPairIntegral Integral_B(Int_B);
PhotoPairIntegral Integral_C(Int_C);

EXPECT_TRUE(Integral_A != Integral_B);
EXPECT_TRUE(Integral_A != Integral_C);

InterpolationDef InterpolDef;
PhotoPairInterpolant PhotoPairInterpolant_A(PhotoPair_A, PhotoAngle_1, InterpolDef);
PhotoPairInterpolant PhotoPairInterpolant_B(PhotoPair_A, PhotoAngle_2, InterpolDef);
PhotoPairInterpolant PhotoPairInterpolant_C(PhotoPair_B, PhotoAngle_1, InterpolDef);
EXPECT_TRUE(PhotoPairInterpolant_A != PhotoPairInterpolant_B);
EXPECT_TRUE(PhotoPairInterpolant_A != PhotoPairInterpolant_C);

PhotoPairInterpolant Interpolant_A(PhotoPairInterpolant_A);
PhotoPairInterpolant Interpolant_B(PhotoPairInterpolant_B);
PhotoPairInterpolant Interpolant_C(PhotoPairInterpolant_C);
EXPECT_TRUE(Integral_A != Interpolant_B);
EXPECT_TRUE(Integral_A != Interpolant_C);

}

TEST(Comparison, Comparison_not_equal_PhotoAngleDistribution)
{
ParticleDef photon  = GammaDef::Get();
ParticleDef mu_def  = MuMinusDef::Get(); //makes no physical sense but programmatically possible
auto medium_1 = std::make_shared<const Water>();
auto medium_2 = std::make_shared<const Ice>();

PhotoAngleNoDeflection PhotoPair_A(photon, medium_1);
PhotoAngleNoDeflection PhotoPair_B(photon, medium_2);
PhotoAngleNoDeflection PhotoPair_C(mu_def, medium_1);
EXPECT_TRUE(PhotoPair_A != PhotoPair_B);
EXPECT_TRUE(PhotoPair_A != PhotoPair_C);

PhotoAngleEGS PhotoPair_D(photon, medium_1);
PhotoAngleEGS PhotoPair_E(photon, medium_2);
PhotoAngleEGS PhotoPair_F(mu_def, medium_1);
EXPECT_TRUE(PhotoPair_D != PhotoPair_E);
EXPECT_TRUE(PhotoPair_D != PhotoPair_F);

PhotoAngleTsaiIntegral PhotoPair_G(photon, medium_1);
PhotoAngleTsaiIntegral PhotoPair_H(photon, medium_2);
PhotoAngleTsaiIntegral PhotoPair_I(mu_def, medium_1);
EXPECT_TRUE(PhotoPair_G != PhotoPair_H);
EXPECT_TRUE(PhotoPair_G != PhotoPair_I);

PhotoAngleDistribution* PhotoAngle_J = new PhotoAngleNoDeflection(photon, medium_1);
PhotoAngleDistribution* PhotoAngle_K = new PhotoAngleTsaiIntegral(photon, medium_1);
EXPECT_TRUE(*PhotoAngle_J != *PhotoAngle_K);
}

TEST(Assignment, Copyconstructor)
{
ParticleDef particle_def = GammaDef::Get();
auto medium = std::make_shared<const Water>();
double multiplier = 1.;

PhotoPairTsai PhotoPair_A(particle_def, medium, multiplier);
PhotoPairTsai PhotoPair_B = PhotoPair_A;
EXPECT_TRUE(PhotoPair_A == PhotoPair_B);

PhotoAngleNoDeflection const PhotoAngle(particle_def, medium);

PhotoPairIntegral Int_A(PhotoPair_A, PhotoAngle);
PhotoPairIntegral Int_B = Int_A;
EXPECT_TRUE(Int_A == Int_B);

InterpolationDef InterpolDef;
PhotoPairInterpolant PhotoPairInterpol_A(PhotoPair_A, PhotoAngle, InterpolDef);
PhotoPairInterpolant PhotoPairInterpol_B = PhotoPairInterpol_A;
EXPECT_TRUE(PhotoPairInterpol_A == PhotoPairInterpol_B);
}

TEST(Assignment, Copyconstructor_PhotoAngleDistribution)
{
ParticleDef particle_def = GammaDef::Get();
auto medium = std::make_shared<const Water>();

PhotoAngleNoDeflection PhotoAngle_A(particle_def, medium);
PhotoAngleNoDeflection PhotoAngle_B = PhotoAngle_A;
EXPECT_TRUE(PhotoAngle_A == PhotoAngle_B);

PhotoAngleTsaiIntegral PhotoAngle_C(particle_def, medium);
PhotoAngleTsaiIntegral PhotoAngle_D = PhotoAngle_C;
EXPECT_TRUE(PhotoAngle_C == PhotoAngle_D);

PhotoAngleEGS PhotoAngle_E(particle_def, medium);
PhotoAngleEGS PhotoAngle_F = PhotoAngle_E;
EXPECT_TRUE(PhotoAngle_E == PhotoAngle_F);
}

TEST(PhotoPair, Test_of_dNdx)
{
std::string filename = "bin/TestFiles/PhotoPair_dNdx.txt";
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

ParticleDef particle_def             = getParticleDef(particleName);
std::shared_ptr<const Medium> medium = CreateMedium(mediumName);

PhotoPairFactory::Definition photopair_def;
photopair_def.multiplier      = multiplier;
photopair_def.parametrization = PhotoPairFactory::Get().GetEnumFromString(parametrization);
photopair_def.photoangle = PhotoPairFactory::PhotoAngle::PhotoAngleNoDeflection;

CrossSection* photopair = PhotoPairFactory::Get().CreatePhotoPair(particle_def, medium, photopair_def);
dNdx_new = photopair->CalculatedNdx(energy);

ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);

delete photopair;
}
}

TEST(PhotoPair, Test_of_dNdx_rnd)
{
std::string filename = "bin/TestFiles/PhotoPair_dNdx_rnd.txt";
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
std::shared_ptr<const Medium> medium = CreateMedium(mediumName);

PhotoPairFactory::Definition photopair_def;
photopair_def.multiplier      = multiplier;
photopair_def.parametrization = PhotoPairFactory::Get().GetEnumFromString(parametrization);
photopair_def.photoangle = PhotoPairFactory::PhotoAngle::PhotoAngleNoDeflection;

CrossSection* photopair = PhotoPairFactory::Get().CreatePhotoPair(particle_def, medium, photopair_def);

dNdx_rnd_new = photopair->CalculatedNdx(energy, rnd);

ASSERT_NEAR(dNdx_rnd_new, dNdx_rnd_stored, 1E-10 * dNdx_rnd_stored);

delete photopair;
}
}

TEST(PhotoPair, Test_of_dNdx_Interpolant)
{
std::string filename = "bin/TestFiles/PhotoPair_dNdx_interpol.txt";
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
std::shared_ptr<const Medium> medium = CreateMedium(mediumName);

PhotoPairFactory::Definition photopair_def;
photopair_def.multiplier      = multiplier;
photopair_def.parametrization = PhotoPairFactory::Get().GetEnumFromString(parametrization);
photopair_def.photoangle = PhotoPairFactory::PhotoAngle::PhotoAngleNoDeflection;

CrossSection* photopair = PhotoPairFactory::Get().CreatePhotoPair(particle_def, medium, photopair_def, InterpolDef);

dNdx_new = photopair->CalculatedNdx(energy);

ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);

delete photopair;
}
}

TEST(PhotoPair, Test_of_dNdxrnd_interpol)
{
std::string filename = "bin/TestFiles/PhotoPair_dNdx_rnd_interpol.txt";
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
std::shared_ptr<const Medium> medium = CreateMedium(mediumName);

PhotoPairFactory::Definition photopair_def;
photopair_def.multiplier      = multiplier;
photopair_def.parametrization = PhotoPairFactory::Get().GetEnumFromString(parametrization);
photopair_def.photoangle = PhotoPairFactory::PhotoAngle::PhotoAngleNoDeflection;

CrossSection* photopair = PhotoPairFactory::Get().CreatePhotoPair(particle_def, medium, photopair_def, InterpolDef);

dNdx_rnd_new = photopair->CalculatedNdx(energy, rnd);

ASSERT_NEAR(dNdx_rnd_new, dNdx_rnd_stored, 1E-10 * dNdx_rnd_stored);

delete photopair;
}
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
