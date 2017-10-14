
#include "gtest/gtest.h"

#include "PROPOSAL/particle/PROPOSALParticle.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

TEST(Comparison , Comparison_equal ) {
    ParticleDef A;
    ParticleDef B;
    EXPECT_TRUE(A==B);

    ParticleDef* C = new ParticleDef();
    ParticleDef* D = new ParticleDef();

    EXPECT_TRUE(*C==*D);
    delete C;
    delete D;

    C = new ParticleDef(ParticleDef::Builder().SetMuMinus().build());
    D = new ParticleDef(ParticleDef::Builder().SetMuMinus().build());

    EXPECT_TRUE(*C==*D);
    delete C;
    delete D;

    C = new ParticleDef(ParticleDef::Builder().SetTauMinus().build());
    D = new ParticleDef(ParticleDef::Builder().SetTauMinus().build());

    EXPECT_TRUE(*C==*D);
    delete C;
    delete D;
}

TEST(Comparison , Comparison_not_equal ) {
    ParticleDef A;
    ParticleDef B;
    EXPECT_TRUE(A==B);

    ParticleDef AA = ParticleDef::Builder().SetMass(100).build();
    EXPECT_TRUE(AA!=B);

    ParticleDef* C = new ParticleDef(ParticleDef::Builder().SetMuMinus().build());
    ParticleDef* D = new ParticleDef(ParticleDef::Builder().SetTauMinus().build());

    EXPECT_TRUE(*C!=*D);
    delete C;
    delete D;
}

TEST(Assignment , Copyconstructor ) {
    ParticleDef A;
    ParticleDef B = A;

    EXPECT_TRUE(A==B);
}

TEST(Assignment , Copyconstructor2 ) {
    ParticleDef A(ParticleDef::Builder().SetMuMinus().build());
    ParticleDef B(A);
    EXPECT_TRUE(A==B);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
