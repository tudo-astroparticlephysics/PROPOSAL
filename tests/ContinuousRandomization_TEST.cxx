
#include "gtest/gtest.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/propagation_utility/ContRandBuilder.h"
#include "PROPOSALTestUtilities/TestFilesHandling.h"

#include "PROPOSAL/crosssection/ParticleDefaultCrossSectionList.h"

using namespace PROPOSAL;

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "MuMinus") {
        return MuMinusDef();
    } else if (name == "TauMinus") {
        return TauMinusDef();
    } else {
        return EMinusDef();
    }
}

TEST(ContinuousRandomization, Constructor)
{
    auto p_def = MuMinusDef();
    auto medium = Ice();
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 1, true);
    auto cross = GetStdCrossSections(p_def, medium, cuts, false);
    auto cont_rand = make_contrand(cross, false);
}

TEST(ContinuousRandomization, FirstMomentum)
{
    // The first momentum of the randomized energy is supposed to be equal to
    // the final_energy when our distribution is neither cut by our upper limit
    // (initial_energy) or our lower limit (mass/low)
    auto p_def = MuMinusDef();
    auto medium = Ice();
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 1, true);
    auto cross = GetStdCrossSections(p_def, medium, cuts, true);
    auto contrand = make_contrand(cross, false);

    RandomGenerator::Get().SetSeed(24601);
    auto energies = std::array<std::pair<double, double>, 5> { { { 5e4, 1e5 },
        { 5e6, 1e7 }, { 5e8, 1e9 }, { 5e10, 1e11 }, { 5e12, 1e13 } } };
    int statistics = 1e4;
    double average, randomized;
    for (auto E : energies) {
        average = 0;
        for (int n = 1; n <= statistics; n++) {
            randomized = contrand->EnergyRandomize(
                E.second, E.first, RandomGenerator::Get().RandomDouble());
            average = average + (randomized - average) / n;
        }
        EXPECT_FALSE(average == E.first);
        EXPECT_NEAR(average, E.first, E.first * 1e-2);
    }
}

TEST(ContinuousRandomization, IdenticalEnergies)
{
    // Expect the randomized energy to be equal to the final energy when the
    // final energy is equal to the initial energy (e.g. the continuous loss is
    // zero)
    auto p_def = MuMinusDef();
    auto medium = Ice();
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 1, true);
    auto cross = GetStdCrossSections(p_def, medium, cuts, true);
    auto contrand = make_contrand(cross, false);

    RandomGenerator::Get().SetSeed(24601);
    double energy = 1e5;
    double energy_randomized = contrand->EnergyRandomize(
        energy, energy, RandomGenerator::Get().RandomDouble());
    EXPECT_DOUBLE_EQ(energy_randomized, energy);
}

TEST(ContinuousRandomization, PhysicalProperties)
{
    // Expecting an increasing variance in randomized_energies when the energy
    // difference between final_energy and an fixed initial energy is increasing
    auto p_def = MuMinusDef();
    auto medium = Ice();
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 1, true);
    auto cross = GetStdCrossSections(p_def, medium, cuts, true);
    auto contrand = make_contrand(cross, false);

    double E_i = 1e12;
    std::array<double, 5> final_energies = { 6e11, 5.5e11, 5e11, 4.5e11, 4e11 };

    RandomGenerator::Get().SetSeed(24601);
    unsigned int statistics = 1e4;
    auto average = std::pair<double, double> { 0., 0. };
    double old_variance = 0;
    for (auto E_f : final_energies) {
        average = { 0., 0. };
        for (unsigned int n = 1; n <= statistics; n++) {
            double sampled = contrand->EnergyRandomize(
                E_i, E_f, RandomGenerator::Get().RandomDouble());
            average = welfords_online_algorithm(
                sampled, n, average.first, average.second);
        }
        EXPECT_TRUE(old_variance < average.second);
        old_variance = average.second;
    }
}

TEST(ContinuousRandomization, Constraints)
{
    // The randomized energy should never be below the particle mass or above
    // the initial_energy
    auto p_def = MuMinusDef();
    auto medium = Ice();
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 1, true);
    auto cross = GetStdCrossSections(p_def, medium, cuts, true);
    auto contrand = make_contrand(cross, false);

    RandomGenerator::Get().SetSeed(24601);
    int statistics = 1e3;
    double randomized;

    for (int n = 1; n < statistics; n++) {
        randomized = contrand->EnergyRandomize(
            1e4, p_def.mass, RandomGenerator::Get().RandomDouble());
        EXPECT_GE(randomized, p_def.mass);
        EXPECT_LT(randomized, 1e4);

        randomized = contrand->EnergyRandomize(
            250, p_def.mass, RandomGenerator::Get().RandomDouble());
        EXPECT_GT(randomized, p_def.mass);
        EXPECT_LT(randomized, 250);

        randomized = contrand->EnergyRandomize(
            110, p_def.mass, RandomGenerator::Get().RandomDouble());
        EXPECT_LT(randomized, 110);
        EXPECT_GE(randomized, p_def.mass);
    }
}

TEST(ContinuousRandomization, compare_integral_interpolant)
{
    auto p_def = MuMinusDef();
    auto medium = Ice();
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 1, true);
    auto cross = GetStdCrossSections(p_def, medium, cuts, true);
    auto contrand_integral = make_contrand(cross, false);
    auto contrand_interpol = make_contrand(cross, true);

    RandomGenerator::Get().SetSeed(24601);
    auto energies = std::array<double, 5> { 1e6, 1e7, 1e8, 1e9, 1e10 };

    for (auto E_i : energies) {
        double rnd = RandomGenerator::Get().RandomDouble();
        double randomized_integral
            = contrand_integral->EnergyRandomize(E_i, 1e5, rnd);
        double randomized_interpol
            = contrand_interpol->EnergyRandomize(E_i, 1e5, rnd);
        EXPECT_NEAR(randomized_integral, randomized_interpol, 1e-3 * E_i);
    }
}

TEST(ContinuousRandomization, Randomize_interpol)
{
    auto in = getTestFiles("continous_randomization.txt");

    char firstLine[256];
    in.getline(firstLine, 256);
    double initial_energy;
    double final_energy;
    double randomized_energy;
    double randomized_energy_new;
    double vcut;
    double ecut;
    std::string mediumName;
    std::string particleName;
    double rnd;
    double energy_old;
    bool first = true;

    std::cout.precision(16);

    RandomGenerator::Get().SetSeed(0);

    while (in.good()) {
        if (first)
            in >> rnd >> particleName >> mediumName >> ecut >> vcut
                >> initial_energy >> final_energy >> randomized_energy;

        first = false;
        energy_old = -1;

        auto p_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);

        // reprouce old behaviour
        if (ecut == -1) {
            ecut = std::numeric_limits<double>::infinity();
        }
        if (vcut == -1) {
            vcut = 1;
        }
        auto cuts = std::make_shared<EnergyCutSettings>(ecut, vcut, true);

        auto cross = GetStdCrossSections(p_def, *medium, cuts, false);
        auto contrand_integral = make_contrand(cross, false);

        while (energy_old < initial_energy) {
            energy_old = initial_energy;
            randomized_energy_new
                = contrand_integral->EnergyRandomize(initial_energy, final_energy, rnd);

            if (initial_energy > 1e6) {
                EXPECT_NEAR(randomized_energy_new, randomized_energy,
                            1e-3 * randomized_energy);
            } // dE2dx for ionization has been underestimated for
                       // previous versions of PROPOSAL

            in >> rnd >> particleName >> mediumName >> ecut >> vcut
                >> initial_energy >> final_energy >> randomized_energy;
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
