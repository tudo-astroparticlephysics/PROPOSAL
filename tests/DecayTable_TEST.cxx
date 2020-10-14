
#include "gtest/gtest.h"

#include "PROPOSAL/decay/DecayTable.h"
#include "PROPOSAL/decay/LeptonicDecayChannel.h"
#include "PROPOSAL/decay/StableChannel.h"
#include "PROPOSAL/decay/TwoBodyPhaseSpace.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/math/RandomGenerator.h"

using namespace PROPOSAL;

ParticleDef mu  = MuMinusDef();
ParticleDef tau = TauMinusDef();
ParticleDef eminus = EMinusDef();
ParticleDef nue = NuMuDef();
ParticleDef nuebar = NuEBarDef();


TEST(Comparison, Comparison_equal)
{
    DecayTable A;
    DecayTable B;
    EXPECT_TRUE(A == B);

    LeptonicDecayChannel x(eminus, nue, nuebar);
    LeptonicDecayChannel y(eminus, nue, nuebar);
    TwoBodyPhaseSpace z(mu, tau);

    A.addChannel(0.5, x);
    B.addChannel(0.5, x);
    EXPECT_TRUE(A == B);

    A.addChannel(0.1, y);
    B.addChannel(0.1, y);
    EXPECT_TRUE(A == B);

    A.addChannel(0.4, z);
    B.addChannel(0.4, z);
    EXPECT_TRUE(A == B);
}

TEST(Comparison, Comparison_not_equal)
{
    DecayTable A;
    DecayTable B;
    EXPECT_TRUE(A == B);

    LeptonicDecayChannel x(eminus, nue, nuebar);
    LeptonicDecayChannel y(eminus, nue, nuebar);
    TwoBodyPhaseSpace z(mu, tau);
    TwoBodyPhaseSpace u(mu, mu);

    A.addChannel(0.5, x);
    EXPECT_TRUE(A != B);

    B.addChannel(0.5, x);
    EXPECT_TRUE(A == B);

    A.addChannel(0.1, z);
    B.addChannel(0.1, u);
    EXPECT_TRUE(A != B);

    A.addChannel(0.4, y);
    B.addChannel(0.4, y);
    EXPECT_TRUE(A != B);
}

TEST(Assignment, Copyconstructor)
{
    DecayTable A;
    DecayTable B = A;

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Copyconstructor2)
{
    DecayTable A;
    DecayTable B(A);

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Operator)
{
    DecayTable A;
    DecayTable B;

    LeptonicDecayChannel x(eminus, nue, nuebar);
    TwoBodyPhaseSpace y(mu, tau);

    A.addChannel(1.0, x);

    EXPECT_TRUE(A != B);

    B = A;

    EXPECT_TRUE(A == B);

    A.addChannel(1.0, y);

    EXPECT_TRUE(A != B);

    B = A;

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Swap)
{
    DecayTable A;
    DecayTable B;
    EXPECT_TRUE(A == B);

    LeptonicDecayChannel x(eminus, nue, nuebar);
    TwoBodyPhaseSpace y(mu, tau);

    A.addChannel(1.0, x);
    A.addChannel(1.0, y);

    DecayTable C = A;
    DecayTable D = B;

    swap(B, A);
    EXPECT_TRUE(C == B);
    EXPECT_TRUE(D == A);
}

TEST(SelectChannel, Muon)
{

    // Leptinic decay channel in muon case
    ParticleDef mu_def = MuMinusDef();
    DecayChannel& dc_muon = mu_def.decay_table.SelectChannel(0.5);

    LeptonicDecayChannel leptonic_channel(eminus, nue, nuebar);

    EXPECT_TRUE(dc_muon == leptonic_channel);
}

TEST(SelectChannel, Electron)
{
    // Leptinic decay channel in electron case
    ParticleDef electron_def = EMinusDef();
    DecayChannel& dc_electron = electron_def.decay_table.SelectChannel(0.5);

    StableChannel stable_channel;

    EXPECT_TRUE(dc_electron == stable_channel);
}

TEST(SelectChannel, Tau)
{
    // tauon decay channels
    ParticleDef tau_def = TauMinusDef();

    int leptonic_count = 0;
    int twobody_count  = 0;
    double random_ch;

    for (int i = 0; i < 1000; ++i)
    {
        random_ch = RandomGenerator::Get().RandomDouble();
        DecayChannel& dc_tau = tau_def.decay_table.SelectChannel(random_ch);

        if (dynamic_cast<LeptonicDecayChannel*>(&dc_tau))
        {
            leptonic_count++;
        } else if (dynamic_cast<TwoBodyPhaseSpace*>(&dc_tau))
        {
            twobody_count++;
        }
    }

    EXPECT_TRUE(leptonic_count > 0);
    EXPECT_TRUE(twobody_count > 0);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
