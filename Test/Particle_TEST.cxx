#include "gtest/gtest.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include <iostream>

TEST(Comparison , Comparison_equal ) {
    PROPOSALParticle A;
    PROPOSALParticle B;
    EXPECT_TRUE(A==B);
    PROPOSALParticle* C = new PROPOSALParticle(ParticleType::TauMinus,1.,1.,1,20,20,1e5,10);
    PROPOSALParticle* D = new PROPOSALParticle(ParticleType::TauMinus,1.,1.,1,20,20,1e5,10);
    C->SetEnergy(1e6);
    D->SetEnergy(1e6);
    EXPECT_TRUE(*C==*D);
    PROPOSALParticle* E = new PROPOSALParticle(ParticleType::MuMinus,0,0,0,0,0,0,0);
    EXPECT_TRUE(A==*E);
    PROPOSALParticle F(1,2,ParticleType::MuMinus,0.,0.,0.,0.,0.,1e5,0.,5,E);
    PROPOSALParticle G(1,2,ParticleType::MuMinus,0.,0.,0.,0.,0.,1e5,0.,5,E);
    EXPECT_TRUE(F==G);
    PROPOSALParticle H(1,2,ParticleType::MuMinus,0.,0.,0.,0.,0.,1e5,0.,5);
    PROPOSALParticle I(1,2,ParticleType::MuMinus,0.,0.,0.,0.,0.,1e5,0.,5);
    EXPECT_TRUE(H==I);
    PROPOSALParticle* J = new PROPOSALParticle(ParticleType::MuMinus);
    PROPOSALParticle* K = new PROPOSALParticle(ParticleType::MuMinus);
    EXPECT_TRUE(*J==*K);


}

TEST(Comparison , Comparison_not_equal ) {
    PROPOSALParticle A;
    PROPOSALParticle B(ParticleType::TauMinus,1.,1.,1,20,20,1e5,10);
    EXPECT_TRUE(A!=B);
    PROPOSALParticle* C = new PROPOSALParticle(ParticleType::TauMinus,1.,1.,1,20,20,1e5,10);
    PROPOSALParticle* D = new PROPOSALParticle(ParticleType::TauMinus,1.,1.,1,20,20,1e5,10);
    D->SetEnergy(1e6);
    EXPECT_TRUE(*C!=*D);

}

TEST(Assignment , Copyconstructor ) {
    PROPOSALParticle A;
    PROPOSALParticle B =A;

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Copyconstructor2 ) {
    PROPOSALParticle A(ParticleType::TauMinus,1.,1.,1,20,20,1e5,10);
    PROPOSALParticle B(A);

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Operator ) {
    PROPOSALParticle A;
    PROPOSALParticle B(ParticleType::TauMinus,1.,1.,1,20,20,1e5,10);

    EXPECT_TRUE(A!=B);

    B=A;

    EXPECT_TRUE(A==B);

    A.SetEnergy(1e12);

    EXPECT_TRUE(A!=B);

    B=A;

    EXPECT_TRUE(A==B);
}

TEST(Assignment , Swap ) {
    PROPOSALParticle A;
    PROPOSALParticle B;
    EXPECT_TRUE(A==B);
    PROPOSALParticle* C = new PROPOSALParticle(ParticleType::TauMinus,1.,1.,1,20,20,1e5,10);
    PROPOSALParticle* D = new PROPOSALParticle(ParticleType::TauMinus,1.,1.,1,20,20,1e5,10);
    C->SetEnergy(1e6);
    D->SetEnergy(1e6);
    EXPECT_TRUE(*C==*D);
    PROPOSALParticle* E = new PROPOSALParticle(ParticleType::MuMinus,0,0,0,0,0,0,0);
    EXPECT_TRUE(A==*E);
    D->swap(A);
    EXPECT_TRUE(*C==A);
    EXPECT_TRUE(*D==B);
    D->swap(*E);
    EXPECT_TRUE(B==*D);

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
