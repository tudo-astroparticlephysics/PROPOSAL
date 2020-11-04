
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

    ParticleState A;
    ParticleState B;
    EXPECT_TRUE(A == B);

    ParticleState* C = new ParticleState();
    C->position = position;
    C->direction = direction;
    ParticleState* D = new ParticleState();
    D->position = position;
    D->direction = direction;

    C->energy = 1e6;
    D->energy = 1e6;
    EXPECT_TRUE(*C == *D);

    position    = Vector3D();
    direction   = Vector3D();
    ParticleState* E = new ParticleState();
    E->position = position;
    E->direction = direction;
    EXPECT_TRUE(A == *E);
}

TEST(Comparison, Comparison_not_equal)
{
    direction.SetSphericalCoordinates(1, 20 * PI / 180., 20 * PI / 180.);
    direction.CalculateCartesianFromSpherical();

    ParticleState A;
    ParticleState B;
    A.SetType(ParticleType::EMinus);
    B.SetType(ParticleType::MuMinus);
    EXPECT_TRUE(A != B);

    ParticleState* C = new ParticleState();
    C->position = position;
    C->direction = direction;
    ParticleState* D = new ParticleState();

    D->position = position;
    D->direction = direction;
    D->energy = 1e6;
    EXPECT_TRUE(*C != *D);
}

TEST(Assignment, Copyconstructor)
{
    ParticleState A;
    ParticleState B = A;

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Copyconstructor2)
{
    direction.SetSphericalCoordinates(1, 20 * PI / 180., 20 * PI / 180.);
    direction.CalculateCartesianFromSpherical();
    ParticleState A;
    A.SetType(ParticleType::TauMinus);
    A.position = position;
    A.direction = direction;
    ParticleState B(A);

    EXPECT_TRUE(A == B);
}

TEST(Deflection, Deflection)
{
    Vector3D direction_A(1., 0., 0.);
    Vector3D direction_B(0., -1./SQRT2, 1./SQRT2);
    Vector3D direction_C(1./3., 2./3., -2./3.);
    Vector3D direction_tmp(1., 0., 0.);
    ParticleState particle;
    particle.SetType(ParticleType::MuMinus);
    double cosangle;

    std::vector<double> cos_phi_list{-1, -0.8, -0.2, 0., 0.2, 0.8, 1.};
    std::vector<double> theta_list{0, PI/4., PI/2., PI, 3./2. * PI, 2. * PI};

    for(auto const& cos_phi: cos_phi_list){
        for(auto const& theta: theta_list){
            particle.direction = direction_A;
            particle.DeflectDirection(cos_phi, theta);
            direction_tmp = particle.direction;
            cosangle = scalar_product(direction_A, direction_tmp);
            ASSERT_NEAR(cos_phi, cosangle, 1e-8);

            particle.direction = direction_B;
            particle.DeflectDirection(cos_phi, theta);
            direction_tmp = particle.direction;
            cosangle = scalar_product(direction_B, direction_tmp);
            ASSERT_NEAR(cos_phi, cosangle, 1e-8);

            particle.direction = direction_C;
            particle.DeflectDirection(cos_phi, theta);
            direction_tmp = particle.direction;
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
