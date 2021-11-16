
#include "gtest/gtest.h"

#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

TEST(Comparison, Comparison_equal)
{
    MuMinusDef A;
    MuMinusDef B;
    EXPECT_TRUE(A == B);

    ParticleDef* C = new MuMinusDef();
    ParticleDef* D = new MuMinusDef();

    EXPECT_TRUE(*C == *D);
    delete C;
    delete D;

    C = new ParticleDef(MuMinusDef());
    D = new ParticleDef(MuMinusDef());

    EXPECT_TRUE(*C == *D);
    delete C;
    delete D;

    C = new ParticleDef(TauMinusDef());
    D = new ParticleDef(TauMinusDef());

    EXPECT_TRUE(*C == *D);
    delete C;
    delete D;
}

TEST(Comparison, Comparison_not_equal)
{
    MuMinusDef A;
    MuMinusDef B;
    EXPECT_TRUE(A == B);

    ParticleDef AA = ParticleDef::Builder().SetMass(100).build();
    EXPECT_TRUE(AA != B);

    ParticleDef* C = new ParticleDef(MuMinusDef());
    ParticleDef* D = new ParticleDef(TauMinusDef());

    EXPECT_TRUE(*C != *D);
    delete C;
    delete D;
}

TEST(Assignment, Copyconstructor)
{
    MuMinusDef A;
    ParticleDef B = A;

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Copyconstructor2)
{
    auto A = MuMinusDef();
    ParticleDef B(A);
    EXPECT_TRUE(A == B);
}

TEST(Type_Particle_Map, GetParticleDefForType)
{
    auto particle1 = EMinusDef();
    EXPECT_EQ(ParticleDef::GetParticleDefForType(particle1.particle_type),
              particle1);
}

TEST(Type_Particle_Map, Registering)
{
    auto particle1 = EMinusDef();
    auto particle1_copy = EMinusDef(); // this should be allowed

    // registering an identical particle is ok
    auto particle2 = ParticleDef::Builder().SetParticleDef(particle1).build();

    // registering a diverging particle with the same id is not allowed
    ASSERT_THROW(ParticleDef::Builder().SetParticleDef(particle1).SetMass(2).build(), std::invalid_argument);

    // registering a diverging particle with an unique id is allowed
    auto particle3 = ParticleDef::Builder().SetParticleDef(particle1).SetMass(2).SetParticleType(100).build();
    EXPECT_EQ(ParticleDef::GetParticleDefForType(particle3.particle_type), particle3);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
