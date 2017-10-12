
// #include <iostream>

#include "gtest/gtest.h"

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

Vector3D position(1.,1.,1.);
Vector3D direction(0.,0.,0.);

TEST(Comparison , Comparison_equal ) {
    direction.SetSphericalCoordinates(1,20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();

    Particle A;
    Particle B;
    EXPECT_TRUE(A==B);

    Particle* C = new Particle(TauMinusDef::Get());
    C->SetPosition(position);
    C->SetDirection(direction);
    Particle* D = new Particle(TauMinusDef::Get());
    D->SetPosition(position);
    D->SetDirection(direction);

    C->SetEnergy(1e6);
    D->SetEnergy(1e6);
    EXPECT_TRUE(*C==*D);

    position = Vector3D();
    direction = Vector3D();
    Particle* E = new Particle(MuMinusDef::Get());
    E->SetPosition(position);
    E->SetDirection(direction);
    EXPECT_TRUE(A==*E);
}

TEST(Comparison , Comparison_not_equal ) {
    direction.SetSphericalCoordinates(1,20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();

    Particle A;
    Particle B(TauMinusDef::Get());
    EXPECT_TRUE(A!=B);

    Particle* C = new Particle(TauMinusDef::Get());
    C->SetPosition(position);
    C->SetDirection(direction);
    Particle* D = new Particle(TauMinusDef::Get());
    D->SetPosition(position);
    D->SetDirection(direction);
    D->SetEnergy(1e6);
    EXPECT_TRUE(*C!=*D);
}

TEST(Assignment , Copyconstructor ) {
    Particle A;
    Particle B = A;

    EXPECT_TRUE(A==B);
}

TEST(Assignment , Copyconstructor2 ) {
    direction.SetSphericalCoordinates(1,20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();
    Particle A(TauMinusDef::Get());
    A.SetPosition(position);
    A.SetDirection(direction);
    Particle B(A);

    EXPECT_TRUE(A==B);
}

// TEST(Assignment , Operator ) {
//     direction.SetSphericalCoordinates(1,20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();
//     Particle A;
//     Particle B(TauMinusDef::Get());
//     B.SetPosition(position);
//     B.SetDirection(direction);
//
//     EXPECT_TRUE(A!=B);
//
//     B=A;
//
//     EXPECT_TRUE(A==B);
//
//     A.SetEnergy(1e12);
//
//     EXPECT_TRUE(A!=B);
//
//     B=A;
//
//     EXPECT_TRUE(A==B);
// }
//
// TEST(Assignment , Swap ) {
//     direction.SetSphericalCoordinates(1,20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();
//     Particle A;
//     Particle B;
//     EXPECT_TRUE(A==B);
//     Particle* C = new Particle(TauMinusDef::Get());
//     C->SetPosition(position);
//     C->SetDirection(direction);
//     Particle* D = new Particle(TauMinusDef::Get());
//     D->SetPosition(position);
//     D->SetDirection(direction);
//
//     C->SetEnergy(1e6);
//     D->SetEnergy(1e6);
//     EXPECT_TRUE(*C==*D);
//
//     position = Vector3D();
//     direction = Vector3D();
//     Particle* E = new Particle(MuMinusDef::Get());
//     E->SetPosition(position);
//     E->SetDirection(direction);
//     EXPECT_TRUE(A==*E);
//
//     D->swap(A);
//     EXPECT_TRUE(*C==A);
//     EXPECT_TRUE(*D==B);
//     D->swap(*E);
//     EXPECT_TRUE(B==*D);
//
// }

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
