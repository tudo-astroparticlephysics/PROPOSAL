
#include "gtest/gtest.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

using namespace PROPOSAL;

Vector3D position(1., 1., 1.);
Vector3D direction(0., 0., 0.);

TEST(Comparison, Comparison_equal)
{
    direction.SetSphericalCoordinates(1, 20 * PI / 180., 20 * PI / 180.);
    direction.CalculateCartesianFromSpherical();

    Particle A;
    Particle B;
    EXPECT_TRUE(A == B);

    Particle* C = new Particle(TauMinusDef::Get());
    C->SetPosition(position);
    C->SetDirection(direction);
    Particle* D = new Particle(TauMinusDef::Get());
    D->SetPosition(position);
    D->SetDirection(direction);

    C->SetEnergy(1e6);
    D->SetEnergy(1e6);
    EXPECT_TRUE(*C == *D);

    position    = Vector3D();
    direction   = Vector3D();
    Particle* E = new Particle(MuMinusDef::Get());
    E->SetPosition(position);
    E->SetDirection(direction);
    EXPECT_TRUE(A == *E);
}

TEST(Comparison, Comparison_not_equal)
{
    direction.SetSphericalCoordinates(1, 20 * PI / 180., 20 * PI / 180.);
    direction.CalculateCartesianFromSpherical();

    Particle A;
    Particle B(TauMinusDef::Get());
    EXPECT_TRUE(A != B);

    Particle* C = new Particle(TauMinusDef::Get());
    C->SetPosition(position);
    C->SetDirection(direction);
    Particle* D = new Particle(TauMinusDef::Get());

    D->SetPosition(position);
    D->SetDirection(direction);
    D->SetEnergy(1e6);
    EXPECT_TRUE(*C != *D);
}

TEST(Assignment, Copyconstructor)
{
    Particle A;
    Particle B = A;

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Copyconstructor2)
{
    direction.SetSphericalCoordinates(1, 20 * PI / 180., 20 * PI / 180.);
    direction.CalculateCartesianFromSpherical();
    Particle A(TauMinusDef::Get());
    A.SetPosition(position);
    A.SetDirection(direction);
    Particle B(A);

    EXPECT_TRUE(A == B);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
