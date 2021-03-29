
#include "gtest/gtest.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/Factories/MupairProductionFactory.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/secondaries/parametrization/mupairproduction/KelnerKokoulinPetrukhinMupairProduction.h"
#include "PROPOSALTestUtilities/TestFilesHandling.h"

using namespace PROPOSAL;

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "MuMinus")
    {
        return MuMinusDef();
    } else if (name == "TauMinus")
    {
        return TauMinusDef();
    } else
    {
        return EMinusDef();
    }
}

const std::string testfile_dir = "tests/TestFiles/";

// TEST(Comparison, Comparison_equal)
// {
// ParticleDef particle_def = MuMinusDef::Get();
// auto medium = std::make_shared<const Water>();
// EnergyCutSettings ecuts;
// double multiplier   = 1.;

// MupairProductionRhoIntegral* MupairInt_A =
//         new MupairKelnerKokoulinPetrukhin(particle_def, medium, ecuts, multiplier, true);
// Parametrization* MupairInt_B = new MupairKelnerKokoulinPetrukhin(particle_def, medium, ecuts, multiplier, true);
// EXPECT_TRUE(*MupairInt_A == *MupairInt_B);

// MupairKelnerKokoulinPetrukhin param_int(particle_def, medium, ecuts, multiplier, true);
// EXPECT_TRUE(param_int == *MupairInt_A);

// MupairIntegral* Int_A        = new MupairIntegral(param_int);
// CrossSectionIntegral* Int_B = new MupairIntegral(param_int);
// EXPECT_TRUE(*Int_A == *Int_B);

// InterpolationDef InterpolDef;
// MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin>* MupairInterpol_A =
//         new MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin>(particle_def, medium, ecuts, multiplier, true, InterpolDef);
// Parametrization* MupairInterpol_B =
//         new MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin>(particle_def, medium, ecuts, multiplier, true, InterpolDef);
// EXPECT_TRUE(*MupairInterpol_A == *MupairInterpol_B);

// MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> param_interpol(particle_def, medium, ecuts, multiplier, true, InterpolDef);
// EXPECT_TRUE(param_interpol == *MupairInterpol_A);

// MupairInterpolant* Interpol_A        = new MupairInterpolant(param_interpol, InterpolDef);
// CrossSectionInterpolant* Interpol_B = new MupairInterpolant(param_interpol, InterpolDef);
// EXPECT_TRUE(*Interpol_A == *Interpol_B);

// delete MupairInt_A;
// delete MupairInt_B;
// delete Int_A;
// delete Int_B;
// delete MupairInterpol_A;
// delete MupairInterpol_B;
// delete Interpol_A;
// delete Interpol_B;
// }

// TEST(Comparison, Comparison_not_equal)
// {
// ParticleDef mu_def  = MuMinusDef::Get();
// ParticleDef tau_def = TauMinusDef::Get();
// auto medium_1 = std::make_shared<const Water>();
// auto medium_2 = std::make_shared<const Ice>();
// EnergyCutSettings ecuts_1(500, -1);
// EnergyCutSettings ecuts_2(-1, 0.05);
// double multiplier_1 = 1.;
// double multiplier_2 = 2.;

// MupairKelnerKokoulinPetrukhin MupairInt_A(mu_def, medium_1, ecuts_1, multiplier_1, true);
// MupairKelnerKokoulinPetrukhin MupairInt_B(tau_def, medium_1, ecuts_1, multiplier_1, true);
// MupairKelnerKokoulinPetrukhin MupairInt_C(mu_def, medium_2, ecuts_1, multiplier_1, true);
// MupairKelnerKokoulinPetrukhin MupairInt_D(mu_def, medium_1, ecuts_2, multiplier_1, true);
// MupairKelnerKokoulinPetrukhin MupairInt_E(mu_def, medium_1, ecuts_1, multiplier_2, true);
// MupairKelnerKokoulinPetrukhin MupairInt_F(mu_def, medium_1, ecuts_1, multiplier_1, false);
// EXPECT_TRUE(MupairInt_A != MupairInt_B);
// EXPECT_TRUE(MupairInt_A != MupairInt_C);
// EXPECT_TRUE(MupairInt_A != MupairInt_D);
// EXPECT_TRUE(MupairInt_A != MupairInt_E);
// EXPECT_TRUE(MupairInt_A != MupairInt_F);

// MupairIntegral Int_A(MupairInt_A);
// MupairIntegral Int_B(MupairInt_B);
// EXPECT_TRUE(Int_A != Int_B);

// InterpolationDef InterpolDef;
// MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_A(mu_def, medium_1, ecuts_1, multiplier_1, true, InterpolDef);
// MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_B(tau_def, medium_1, ecuts_1, multiplier_1, true, InterpolDef);
// MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_C(mu_def, medium_2, ecuts_1, multiplier_1, true, InterpolDef);
// MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_D(mu_def, medium_1, ecuts_2, multiplier_1, true, InterpolDef);
// MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_E(mu_def, medium_1, ecuts_1, multiplier_2, true, InterpolDef);
// MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_F(mu_def, medium_1, ecuts_1, multiplier_1, false, InterpolDef);
// EXPECT_TRUE(MupairInterpol_A != MupairInterpol_B);
// EXPECT_TRUE(MupairInterpol_A != MupairInterpol_C);
// EXPECT_TRUE(MupairInterpol_A != MupairInterpol_D);
// EXPECT_TRUE(MupairInterpol_A != MupairInterpol_E);
// EXPECT_TRUE(MupairInterpol_A != MupairInterpol_F);

// MupairInterpolant Interpol_A(MupairInterpol_A, InterpolDef);
// MupairInterpolant Interpol_B(MupairInterpol_B, InterpolDef);
// EXPECT_TRUE(Interpol_A != Interpol_B);
// }

// TEST(Assignment, Copyconstructor)
// {
// ParticleDef particle_def = MuMinusDef::Get();
// auto medium = std::make_shared<const Water>();
// EnergyCutSettings ecuts;
// double multiplier = 1.;

// MupairKelnerKokoulinPetrukhin MupairInt_A(particle_def, medium, ecuts, true, multiplier);
// MupairKelnerKokoulinPetrukhin MupairInt_B = MupairInt_A;
// EXPECT_TRUE(MupairInt_A == MupairInt_B);

// MupairIntegral Int_A(MupairInt_A);
// MupairIntegral Int_B = Int_A;
// EXPECT_TRUE(Int_A == Int_B);

// InterpolationDef InterpolDef;
// MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_A(particle_def, medium, ecuts, multiplier, true, InterpolDef);
// MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_B = MupairInterpol_A;
// EXPECT_TRUE(MupairInterpol_A == MupairInterpol_B);

// MupairInterpolant Interpol_A(MupairInterpol_A, InterpolDef);
// MupairInterpolant Interpol_B = Interpol_A;
// EXPECT_TRUE(Interpol_A == Interpol_B);
// }

// TEST(Assignment, Copyconstructor2)
// {
// ParticleDef particle_def = MuMinusDef::Get();
// auto medium = std::make_shared<const Water>();
// EnergyCutSettings ecuts;
// double multiplier = 1.;

// MupairKelnerKokoulinPetrukhin MupairInt_A(particle_def, medium, ecuts, true, multiplier);
// MupairKelnerKokoulinPetrukhin MupairInt_B(MupairInt_A);
// EXPECT_TRUE(MupairInt_A == MupairInt_B);

// MupairIntegral Int_A(MupairInt_A);
// MupairIntegral Int_B(Int_A);
// EXPECT_TRUE(Int_A == Int_B);

// InterpolationDef InterpolDef;
// MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_A(particle_def, medium, ecuts, multiplier, true, InterpolDef);
// MupairProductionRhoInterpolant<MupairKelnerKokoulinPetrukhin> MupairInterpol_B(MupairInterpol_A);
// EXPECT_TRUE(MupairInterpol_A == MupairInterpol_B);

// MupairInterpolant Interpol_A(MupairInterpol_A, InterpolDef);
// MupairInterpolant Interpol_B(Interpol_A);
// EXPECT_TRUE(Interpol_A == Interpol_B);
// }

// in polymorphism an assignment and swap operator doesn't make sense

TEST(Mupairproduction, Test_of_dEdx)
{
    auto in = getTestFiles("Mupair_dEdx.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    double energy;
    std::string parametrization;
    double dEdx_stored;
    double dEdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> parametrization >> dEdx_stored)
    {
        parametrization.erase(0,6);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        if (particleName != "MuMinus")
            continue; // parametrization at the moment only optimized for ingoing muons

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_mupairproduction(particle_def, *medium, ecuts, false,
                                           config);

        dEdx_new = cross->CalculatedEdx(energy) * medium->GetMassDensity();
        if (energy <= 10000)
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-1 * dEdx_stored); // there have been changes in the cross section parametrization
        else
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);
    }
}

TEST(Mupairproduction, Test_of_dNdx)
{
    auto in = getTestFiles("Mupair_dNdx.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> parametrization >> dNdx_stored)
    {
        parametrization.erase(0,6);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        if (particleName != "MuMinus")
            continue; // parametrization at the moment only optimized for ingoing muons

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_mupairproduction(particle_def, *medium, ecuts, false,
                                           config);

        dNdx_new = cross->CalculatedNdx(energy) * medium->GetMassDensity();

        if (energy <= 10000)
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-2 * dNdx_stored); // there have been changes in the cross section parametrization
        else
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
    }
}

TEST(Mupairproduction, Test_Stochastic_Loss)
{
    auto in = getTestFiles("Mupair_e.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    double energy;
    std::string parametrization;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> parametrization >> rnd1 >> rnd2 >> stochastic_loss_stored)
    {
        parametrization.erase(0,6);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        if (particleName != "MuMinus")
            continue; // parametrization at the moment only optimized for ingoing muons

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_mupairproduction(particle_def, *medium, ecuts, false,
                                           config);

        auto dNdx_full = cross->CalculatedNdx(energy);
        auto components = medium->GetComponents();
        double sum = 0;

        for (auto comp : components)
        {
            double dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());
            sum += dNdx_for_comp;
            if (sum > dNdx_full * rnd2) {
                double rate_new = dNdx_for_comp * rnd1;
                if (ecut == INF and vcut == 1) {
                    #ifndef NDEBUG
                    EXPECT_DEATH(cross->CalculateStochasticLoss(comp.GetHash(), energy, rate_new), "");
                    #endif
                } else {
                    auto v = cross->CalculateStochasticLoss(comp.GetHash(), energy, rate_new);
                    if (energy <= 10000)
                        EXPECT_NEAR(energy * v, stochastic_loss_stored, 1E-2 * stochastic_loss_stored); // there have been changes in the cross section parametrization
                    else
                        EXPECT_NEAR(energy * v, stochastic_loss_stored, 1E-3 * stochastic_loss_stored); // there have been changes in the cross section parametrization

                    // cross_check
                    auto rate_rnd = cross->CalculateCumulativeCrosssection(energy, comp.GetHash(), v);
                    EXPECT_NEAR(rate_rnd / dNdx_for_comp, rnd1, 1e-3);
                    break;
                }
            }
        }
    }
}

TEST(Mupairproduction, Test_Calculate_Rho)
{
auto in = getTestFiles("Mupair_rho.txt");

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

double rho;


std::cout.precision(16);
RandomGenerator::Get().SetSeed(0);

while (in >> particleName >> mediumName >> ecut >> vcut >> v >> multiplier >> energy >> parametrization >> rnd1 >> rnd2 >> E1_stored >> E2_stored)
{
    if (vcut == -1)
        vcut = 1;
    if (ecut == -1)
        ecut = INF;

if (particleName != "MuMinus")
    continue; // parametrization at the moment only optimized for ingoing muons

ParticleDef particle_def = getParticleDef(particleName);
std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);
auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, false);

auto fac = secondaries::KelnerKokoulinPetrukhinMupairProduction(particle_def, *medium);
rho = fac.CalculateRho(energy, v, medium->GetComponents().front(), rnd1, rnd2);

E1_new = 0.5 * v * energy * (1 + rho);
E2_new = 0.5 * v * energy * (1 - rho);


if (energy <= 1e5) {
    EXPECT_NEAR(E1_new, E1_stored, 5E-2 * E1_stored);
    EXPECT_NEAR(E2_new, E2_stored, 5E-2 * E2_stored);
} else {
    EXPECT_NEAR(E1_new, E1_stored, 1E-3 * E1_stored);
    EXPECT_NEAR(E2_new, E2_stored, 1E-3 * E2_stored);
}
}
}

TEST(Mupairproduction, Test_of_dEdx_Interpolant)
{
    auto in = getTestFiles("Mupair_dEdx.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    double energy;
    std::string parametrization;
    double dEdx_stored;
    double dEdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> parametrization >> dEdx_stored)
    {
        parametrization.erase(0,6);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        if (particleName != "MuMinus")
            continue; // parametrization at the moment only optimized for ingoing muons

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_mupairproduction(particle_def, *medium, ecuts, true,
                                           config);
        auto cross_integral = make_mupairproduction(particle_def, *medium, ecuts,
                                                    false, config);

        dEdx_new = cross->CalculatedEdx(energy) * medium->GetMassDensity();
        auto dEdx_integral = cross_integral->CalculatedEdx(energy) * medium->GetMassDensity();

        if (energy * vcut == ecut) {
            EXPECT_NEAR(dEdx_new, dEdx_stored, 2e-1 * dEdx_stored);
            EXPECT_NEAR(dEdx_new, dEdx_integral, 2e-1 * dEdx_integral);
        } else if (energy <= 10000) {
            EXPECT_NEAR(dEdx_new, dEdx_stored, 5e-2 * dEdx_stored);
            EXPECT_NEAR(dEdx_new, dEdx_integral, 5e-2 * dEdx_integral);
        } else {
            EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);
            EXPECT_NEAR(dEdx_new, dEdx_integral, 1e-3 * dEdx_integral);
        }

        // cross check
    }
}

TEST(Mupairproduction, Test_of_dNdx_Interpolant)
{
    auto in = getTestFiles("Mupair_dNdx.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> parametrization >> dNdx_stored)
    {
        parametrization.erase(0,6);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        if (particleName != "MuMinus")
            continue; // parametrization at the moment only optimized for ingoing muons

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_mupairproduction(particle_def, *medium, ecuts, true,
                                           config);
        auto cross_integral = make_mupairproduction(particle_def, *medium, ecuts,
                                                    false, config);
        dNdx_new = cross->CalculatedNdx(energy) * medium->GetMassDensity();
        auto dNdx_integral = cross_integral->CalculatedNdx(energy) * medium->GetMassDensity();

        if (energy * vcut == ecut) {
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-1 * dNdx_stored); // kink in integrand
            EXPECT_NEAR(dNdx_new, dNdx_integral, 1e-1 * dNdx_integral);
        } else if (energy <= 10000) {
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-2 * dNdx_stored); // parametrization changed for small energies
            EXPECT_NEAR(dNdx_new, dNdx_integral, 1e-3 * dNdx_integral);
        } else {
            EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
            EXPECT_NEAR(dNdx_new, dNdx_integral, 1e-3 * dNdx_integral);
        }

    }
}

TEST(Mupairproduction, Test_of_e_interpol)
{
    auto in = getTestFiles("Mupair_e.txt");

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    double energy;
    std::string parametrization;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;

    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >>  energy >> parametrization >> rnd1 >> rnd2 >> stochastic_loss_stored)
    {
        parametrization.erase(0,6);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        if (particleName != "MuMinus")
            continue; // parametrization at the moment only optimized for ingoing muons

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        nlohmann::json config;
        config["parametrization"] = parametrization;

        auto cross = make_mupairproduction(particle_def, *medium, ecuts, true,
                                           config);

        auto dNdx_full = cross->CalculatedNdx(energy);
        auto components = medium->GetComponents();
        double sum = 0;

        for (auto comp : components)
        {
            double dNdx_for_comp = cross->CalculatedNdx(energy, comp.GetHash());
            sum += dNdx_for_comp;
            if (sum > dNdx_full * rnd2) {
                double rate_new = dNdx_for_comp * rnd1;
                if (ecut == INF and vcut == 1 ) {
                    #ifndef NDEBUG
                    EXPECT_DEATH(cross->CalculateStochasticLoss(comp.GetHash(), energy, rate_new), "");
                    #endif
                } else {
                    auto v = cross->CalculateStochasticLoss(comp.GetHash(), energy, rate_new);
                    if (energy <= 10000)
                        EXPECT_NEAR(v * energy, stochastic_loss_stored, 1E-2 * stochastic_loss_stored); // parametrization changed for small energies
                    else
                        EXPECT_NEAR(v * energy, stochastic_loss_stored, 5E-3 * stochastic_loss_stored);

                    // cross check
                    auto rate_rnd = cross->CalculateCumulativeCrosssection(energy, comp.GetHash(), v);
                    EXPECT_NEAR(rate_rnd/dNdx_for_comp, rnd1, 1e-4);
                    break;
                }
            }
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
