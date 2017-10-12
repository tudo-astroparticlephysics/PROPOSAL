
#include "gtest/gtest.h"

#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

TEST(Comparison , Comparison_equal ) {
    ParticleDef A;
    ParticleDef B;
    EXPECT_TRUE(A==B);

    A.mass = 100;
    B.mass = 100;
    EXPECT_TRUE(A==B);

    ParticleDef* C = new ParticleDef();
    ParticleDef* D = new ParticleDef();

    EXPECT_TRUE(*C==*D);
    delete C;
    delete D;

    C = new ParticleDef(MuMinusDef::Get());
    D = new ParticleDef(MuMinusDef::Get());

    EXPECT_TRUE(*C==*D);
    delete C;
    delete D;

    C = new ParticleDef(TauMinusDef::Get());
    D = new ParticleDef(TauMinusDef::Get());

    EXPECT_TRUE(*C==*D);
    delete C;
    delete D;
}

TEST(Comparison , Comparison_not_equal ) {
    ParticleDef A;
    ParticleDef B;
    EXPECT_TRUE(A==B);

    A.mass = 100;
    EXPECT_TRUE(A!=B);

    ParticleDef* C = new ParticleDef();
    ParticleDef* D = new ParticleDef();
    C->lifetime = 2.0;

    EXPECT_TRUE(*C!=*D);
    delete C;
    delete D;

    C = new ParticleDef(MuMinusDef::Get());
    D = new ParticleDef(TauMinusDef::Get());

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
    ParticleDef A(TauMinusDef::Get());
    ParticleDef B(A);
    EXPECT_TRUE(A==B);
}

TEST(Assignment , Operator ) {
    ParticleDef A;
    ParticleDef B(TauMinusDef::Get());
    std::cout << A << std::endl;
    std::cout << B << std::endl;


    EXPECT_TRUE(A!=B);

    B=A;

    EXPECT_TRUE(A==B);

    A.charge = 2.0;

    EXPECT_TRUE(A!=B);

    B=A;

    EXPECT_TRUE(A==B);
}

TEST(Assignment , Swap ) {
    ParticleDef A(MuMinusDef::Get());
    ParticleDef B(MuMinusDef::Get());
    EXPECT_TRUE(A==B);
    ParticleDef* C = new ParticleDef(TauMinusDef::Get());
    ParticleDef* D = new ParticleDef(TauMinusDef::Get());

    C->lifetime = 2.0;
    D->lifetime = 2.0;
    EXPECT_TRUE(*C==*D);

    ParticleDef* E = new ParticleDef(MuMinusDef::Get());
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
