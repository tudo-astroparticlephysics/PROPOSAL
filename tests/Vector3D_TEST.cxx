
#include <cmath>
#include <iostream>
#include <PROPOSAL/math/Cartesian3D.h>
#include <PROPOSAL/math/Spherical3D.h>

#include "gtest/gtest.h"

#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

TEST(Comparison, Comparison_equal)
{
    Cartesian3D A;
    Cartesian3D B;
    EXPECT_TRUE(A == B);
    Vector3D* C = new Cartesian3D(1., 2., 3.);
    Vector3D* D = new Cartesian3D(1., 2., 3.);
    EXPECT_TRUE(*C == *D);
    D->SetCoordinates({0., 0., 0.});
    EXPECT_TRUE(A == *D);
    B.SetCoordinates({0.1, 0.2, 0.3});
    D->SetCoordinates({0.1, 0.2, 0.3});
    EXPECT_TRUE(B == *D);
}

TEST(Comparison, Comparison_not_equal)
{
    Cartesian3D A;
    Cartesian3D B;
    B.SetCoordinates({1., 2., 3.});
    EXPECT_TRUE(A != B);
    Cartesian3D C;
    C.SetCoordinates({0.1, 0.2, 0.3});
    EXPECT_TRUE(A != C);
}

TEST(Comparison, Comparison_equal_spherical)
{
    Spherical3D A;
    Spherical3D B;
    EXPECT_TRUE(A == B);
    Vector3D* C = new Spherical3D(1., 2., 3.);
    Vector3D* D = new Spherical3D(1., 2., 3.);
    EXPECT_TRUE(*C == *D);
    D->SetCoordinates({0., 0., 0.});
    EXPECT_TRUE(A == *D);
    B.SetCoordinates({0.1, 0.2, 0.3});
    D->SetCoordinates({0.1, 0.2, 0.3});
    EXPECT_TRUE(B == *D);
}

TEST(Comparison, Comparison_not_equal_spherical)
{
    Spherical3D A;
    Spherical3D B;
    B.SetCoordinates({1., 2., 3.});
    EXPECT_TRUE(A != B);
    Spherical3D C;
    C.SetCoordinates({0.1, 0.2, 0.3});
    EXPECT_TRUE(A != C);
}

TEST(Assignment, Copyconstructor)
{
    Cartesian3D A(1, 2, 3);
    Cartesian3D B(A);
    EXPECT_TRUE(A == B);
    Cartesian3D C(2, 3, 4);
    B = C;
    EXPECT_TRUE(C == B);
}

TEST(Assignment, Operator)
{
    Cartesian3D A;
    Cartesian3D B;
    A.SetCoordinates({1, 2, 3});
    EXPECT_TRUE(A != B);
    B = A;
    EXPECT_TRUE(A == B);
}

TEST(Assignment, Copyconstructor_spherical)
{
    Spherical3D A(1, 2, 3);
    Spherical3D B(A);
    EXPECT_TRUE(A == B);
    Spherical3D C(2, 3, 4);
    B = C;
    EXPECT_TRUE(C == B);
}

TEST(Assignment, Operator_spherical)
{
    Spherical3D A;
    Spherical3D B;
    A.SetCoordinates({1, 2, 3});
    EXPECT_TRUE(A != B);
    B = A;
    EXPECT_TRUE(A == B);
}

TEST(Addition, Operator)
{
    Cartesian3D A;
    Cartesian3D B;
    Cartesian3D C;
    Cartesian3D  D;
    A.SetCoordinates({1, 2, 3});
    B.SetCoordinates({2, 4, 6});
    C.SetCoordinates({3, 6, 9});
    D = A + B;
    EXPECT_TRUE(C == D);
    D.SetCoordinates({0, 0, 0});
    EXPECT_TRUE(D != C);
    D = B + A;
    EXPECT_TRUE(D == C);
}

TEST(Subtraction, Operator)
{
    Cartesian3D A;
    Cartesian3D B;
    Cartesian3D C;
    Cartesian3D D;
    Cartesian3D E;
    A.SetCoordinates({1, 2, 3});
    B.SetCoordinates({2, 4, 6});
    C.SetCoordinates({3, 6, 9});
    E.SetCoordinates({-1, -2, -3});
    D = -A;
    EXPECT_TRUE(D == E);
    D = C - B;
    EXPECT_TRUE(D == A);
    D = B - C;
    EXPECT_TRUE(D == E);
}

TEST(Scaling, Operator)
{
    Cartesian3D A;
    Cartesian3D B;
    Cartesian3D C;
    double factor1 = 2.;
    A.SetCoordinates({1, 2, 3});
    B.SetCoordinates({2, 4, 6});
    C = A * factor1;
    EXPECT_TRUE(C == B);
    C.SetCoordinates({0, 0, 0});
    EXPECT_TRUE(C != B);
    C = factor1 * A;
    EXPECT_TRUE(C == B);
}

TEST(ScalarProduct, Operator)
{
    Cartesian3D A;
    Cartesian3D B;
    Cartesian3D C;
    double factor1 = 5;
    double result  = 0;
    A.SetCoordinates({1, 1, 1});
    B.SetCoordinates({2, 1, 2});
    EXPECT_TRUE(result != factor1);
    result = A*B;
    EXPECT_TRUE(result == factor1);
    result = 0.;
    result = B*A;
    EXPECT_TRUE(result == factor1);
}

TEST(VectorProduct, Operator)
{
    Cartesian3D A;
    Cartesian3D B;
    Cartesian3D C;
    Cartesian3D D;
    Cartesian3D E;
    A.SetCoordinates({1, 0, 0});
    B.SetCoordinates({0, 1, 0});
    C.SetCoordinates({0, 0, 1});
    D = vector_product(A, B);
    EXPECT_TRUE(D == C);
    D = -vector_product(B, A);
    EXPECT_TRUE(D == C);
    D = -vector_product(A, C);
    EXPECT_TRUE(D == B);
    D = vector_product(B, C);
    EXPECT_TRUE(D == A);
}

TEST(Magnitude, Operator)
{
    Cartesian3D A;
    Cartesian3D B;
    double sum = 3;
    double result;
    A.SetCoordinates({1, 2, 2});
    B.SetCoordinates({1, 1, 1});
    result = A.magnitude();
    EXPECT_TRUE(result == sum);
    result = B.magnitude();
    EXPECT_TRUE(result != sum);
}

TEST(Magnitude, Operator_spherical)
{
    Spherical3D A;
    Spherical3D B;
    double sum = 1;
    double result;
    A.SetCoordinates({1, 2, 2});
    B.SetCoordinates({3, 1, 1});
    result = A.magnitude();
    EXPECT_TRUE(result == sum);
    result = B.magnitude();
    EXPECT_TRUE(result != sum);
}

TEST(Normalise, Operator)
{
    Cartesian3D A;
    Cartesian3D B;
    Cartesian3D C;
    A.SetCoordinates({3, 0, 0});
    B.SetCoordinates({1, 1, 1});
    C.SetCoordinates({1, 0, 0});
    EXPECT_TRUE(A != C);
    EXPECT_TRUE(B != C);
    A.normalize();
    EXPECT_TRUE(A == C);
    B.normalize();
    EXPECT_TRUE(B != C);
}

TEST(Normalise, Operator_spherical)
{
    Spherical3D A;
    Spherical3D B;
    Spherical3D C;
    A.SetCoordinates({3, 0, 0});
    B.SetCoordinates({1, 1, 1});
    C.SetCoordinates({1, 0, 0});
    EXPECT_TRUE(A != C);
    EXPECT_TRUE(B != C);
    A.normalize();
    EXPECT_TRUE(A == C);
    B.normalize();
    EXPECT_TRUE(B != C);
}

TEST(CalculateSphericalCoordinates, Conversion)
{
    Cartesian3D A(1, 2, 2);
    Spherical3D B(A.GetSphericalCoordinates());
    Spherical3D C(3, std::atan2(2., 1.), std::acos(2. / 3.));
    EXPECT_FALSE(A == B);
    EXPECT_FALSE(A == C);
    EXPECT_TRUE(B == C);
}

TEST(CalculateCartesianFromSpherical, Conversion)
{
    Cartesian3D A(1, 2, 2);
    Spherical3D B(3, std::atan2(2., 1.), std::acos(2. / 3.));
    Cartesian3D C(B);
    double epsilon = std::numeric_limits<double>::epsilon();
    double error_factor = 2.;
    bool test_x = std::abs(A.GetX() - C.GetX()) < std::abs(std::min(A.GetX(), C.GetX())) * epsilon * error_factor;
    bool test_y = std::abs(A.GetY() - C.GetY()) < std::abs(std::min(A.GetY(), C.GetY())) * epsilon * error_factor;
    bool test_z = std::abs(A.GetZ() - C.GetZ()) < std::abs(std::min(A.GetZ(), C.GetZ())) * epsilon * error_factor;
    EXPECT_TRUE(A == C || (test_x && test_y && test_z));
}

TEST(Deflection, Deflection)
{
    Cartesian3D direction_A(1., 0., 0.);
    Cartesian3D direction_B(0., -1./SQRT2, 1./SQRT2);
    Cartesian3D direction_C(1./3., 2./3., -2./3.);
    Cartesian3D direction_tmp(1., 0., 0.);
    double cosangle;

    std::vector<double> cos_phi_list{-1, -0.8, -0.2, 0., 0.2, 0.8, 1.};
    std::vector<double> theta_list{0, PI/4., PI/2., PI, 3./2. * PI, 2. * PI};

    for(auto const& cos_phi: cos_phi_list){
        for(auto const& theta: theta_list){
            direction_tmp = direction_A;
            direction_tmp.deflect(cos_phi, theta);
            cosangle = direction_A*direction_tmp;
            EXPECT_NEAR(cos_phi, cosangle, 1e-8);

            direction_tmp = direction_B;
            direction_tmp.deflect(cos_phi, theta);
            cosangle = direction_B*direction_tmp;
            EXPECT_NEAR(cos_phi, cosangle, 1e-8);

            direction_tmp = direction_C;
            direction_tmp.deflect(cos_phi, theta);
            cosangle = direction_C*direction_tmp;
            EXPECT_NEAR(cos_phi, cosangle, 1e-8);
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
