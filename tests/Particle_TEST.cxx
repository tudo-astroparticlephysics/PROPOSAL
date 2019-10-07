
#include "gtest/gtest.h"
#include <cmath>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/math/Vector3D.h"

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

TEST(Deflection, Deflection)
{
    Vector3D direction_A(1., 0., 0.);
    Vector3D direction_B(0., -1./SQRT2, 1./SQRT2);
    Vector3D direction_C(1./3., 2./3., -2./3.);
    Vector3D direction_tmp(1., 0., 0.);
    Particle particle(MuMinusDef::Get());
    double cosangle;

    std::vector<double> cos_phi_list{-1, -0.8, -0.2, 0., 0.2, 0.8, 1.};
    std::vector<double> theta_list{0, PI/4., PI/2., PI, 3./2. * PI, 2. * PI};

    for(auto const& cos_phi: cos_phi_list){
        for(auto const& theta: theta_list){
            particle.SetDirection(direction_A);
            particle.DeflectDirection(cos_phi, theta);
            direction_tmp = particle.GetDirection();
            cosangle = scalar_product(direction_A, direction_tmp);
            ASSERT_NEAR(cos_phi, cosangle, 1e-8);

            particle.SetDirection(direction_B);
            particle.DeflectDirection(cos_phi, theta);
            direction_tmp = particle.GetDirection();
            cosangle = scalar_product(direction_B, direction_tmp);
            ASSERT_NEAR(cos_phi, cosangle, 1e-8);

            particle.SetDirection(direction_C);
            particle.DeflectDirection(cos_phi, theta);
            direction_tmp = particle.GetDirection();
            cosangle = scalar_product(direction_C, direction_tmp);
            ASSERT_NEAR(cos_phi, cosangle, 1e-8);
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
