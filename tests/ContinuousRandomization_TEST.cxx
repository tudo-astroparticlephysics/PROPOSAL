
#include "gtest/gtest.h"

#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/propagation_utility/ContinuousRandomizer.h"
#include <fstream>

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

class Test_Utilities : public ::testing::Test
{
protected:
    Test_Utilities() {}
    virtual ~Test_Utilities() {}

    static Utility a;
    static Utility b;
};

Utility Test_Utilities::a(MuMinusDef::Get(), Ice(), EnergyCutSettings(), Utility::Definition(), InterpolationDef());
Utility Test_Utilities::b(TauMinusDef::Get(), Ice(), EnergyCutSettings(), Utility::Definition(), InterpolationDef());

TEST_F(Test_Utilities, Comparison_equal)
{
    ContinuousRandomizer A(a);
    ContinuousRandomizer B(a);

    EXPECT_TRUE(A == B);

    ContinuousRandomizer* C = new ContinuousRandomizer(a);
    ContinuousRandomizer* D = new ContinuousRandomizer(a);

    EXPECT_TRUE(*C == *D);

    delete C;
    delete D;

    A.Randomize(0, 1, 0.5);
    B.Randomize(0, 1, 0.5);

    EXPECT_TRUE(A == B);
}

TEST_F(Test_Utilities, Comparison_not_equal)
{
    ContinuousRandomizer A(a);
    ContinuousRandomizer B(b);

    // Only check different particle. If test does not fail, it should be good, because interanlly the utilitys are
    // compared. So the Utility test must work fine.
    EXPECT_TRUE(A != B);
}

TEST_F(Test_Utilities, Copyconstructor)
{
    ContinuousRandomizer A(a);
    ContinuousRandomizer B = A;

    EXPECT_TRUE(A == B);
}

TEST_F(Test_Utilities, Copyconstructor2)
{
    ContinuousRandomizer A(a);
    ContinuousRandomizer B(A);

    EXPECT_TRUE(A == B);

    Utility utility_c(a);

    EXPECT_TRUE(a == utility_c);

    ContinuousRandomizer C(a);
    ContinuousRandomizer D(utility_c, C);

    EXPECT_TRUE(A == B);
}

TEST(ContinuousRandomization, Randomize_interpol)
{
    std::ifstream in;
    std::string filename = "bin/TestFiles/continous_randomization.txt";
    in.open(filename.c_str());

    if (!in.good())
    {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

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

    while (in.good())
    {
        if (first)
            in >> rnd >> particleName >> mediumName >> ecut >> vcut >> initial_energy >> final_energy >>
                randomized_energy;

        first      = false;
        energy_old = -1;

        Medium* medium = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings cut_settings(ecut, vcut);
        ParticleDef particle_def = getParticleDef(particleName);

        Utility utility(particle_def, *medium, cut_settings, Utility::Definition(), InterpolationDef());
        ContinuousRandomizer cont(utility, InterpolationDef());

        while (energy_old < initial_energy)
        {
            energy_old = initial_energy;

            randomized_energy_new = cont.Randomize(initial_energy, final_energy, rnd);

            ASSERT_NEAR(randomized_energy_new, randomized_energy, 1e-1 * randomized_energy);

            in >> rnd >> particleName >> mediumName >> ecut >> vcut >> initial_energy >> final_energy >>
                randomized_energy;
        }

        delete medium;
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
