
#include "gtest/gtest.h"

#include "PROPOSAL/decay/DecayTable.h"
#include "PROPOSAL/decay/LeptonicDecayChannel.h"
#include "PROPOSAL/decay/StableChannel.h"
#include "PROPOSAL/decay/TwoBodyPhaseSpace.h"
#include "PROPOSAL/particle/Particle.h"

using namespace PROPOSAL;

ParticleDef mu  = MuMinusDef::Get();
ParticleDef tau = TauMinusDef::Get();

TEST(Comparison, Comparison_equal)
{
    DecayTable A;
    DecayTable B;
    EXPECT_TRUE(A == B);

    LeptonicDecayChannel x(EMinusDef::Get(), NuEDef::Get(), NuEBarDef::Get());
    LeptonicDecayChannel y(EMinusDef::Get(), NuEDef::Get(), NuEBarDef::Get());
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

    LeptonicDecayChannel x(EMinusDef::Get(), NuEDef::Get(), NuEBarDef::Get());
    LeptonicDecayChannel y(EMinusDef::Get(), NuEDef::Get(), NuEBarDef::Get());
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

    LeptonicDecayChannel x(EMinusDef::Get(), NuEDef::Get(), NuEBarDef::Get());
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

    LeptonicDecayChannel x(EMinusDef::Get(), NuEDef::Get(), NuEBarDef::Get());
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
    Particle muon;
    DecayChannel& dc_muon = muon.GetDecayTable().SelectChannel();

    LeptonicDecayChannel leptonic_channel(EMinusDef::Get(), NuMuDef::Get(), NuEBarDef::Get());

    EXPECT_TRUE(dc_muon == leptonic_channel);
}

TEST(SelectChannel, Electron)
{
    // Leptinic decay channel in electron case
    Particle electron(EMinusDef::Get());
    DecayChannel& dc_electron = electron.GetDecayTable().SelectChannel();

    StableChannel stable_channel;

    EXPECT_TRUE(dc_electron == stable_channel);
}

TEST(SelectChannel, Tau)
{
    // tauon decay channels
    Particle tau(TauMinusDef::Get());

    int leptonic_count = 0;
    int twobody_count  = 0;

    for (int i = 0; i < 1000; ++i)
    {
        DecayChannel& dc_tau = tau.GetDecayTable().SelectChannel();

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
