
#include <cmath>
#include <iostream>

#include "gtest/gtest.h"

#include "PROPOSAL/math/Vector3D.h"

using namespace PROPOSAL;

TEST(Comparison, Comparison_equal)
{
    Vector3D A;
    Vector3D B;
    EXPECT_TRUE(A == B);
    Vector3D* C = new Vector3D(1., 2., 3.);
    Vector3D* D = new Vector3D(1., 2., 3.);
    EXPECT_TRUE(*C == *D);
    D->SetCartesianCoordinates(0., 0., 0.);
    EXPECT_TRUE(A == *D);
    B.SetSphericalCoordinates(0.1, 0.2, 0.3);
    D->SetSphericalCoordinates(0.1, 0.2, 0.3);
    EXPECT_TRUE(B == *D);
}

TEST(Comparison, Comparison_not_equal)
{
    Vector3D A;
    Vector3D B;
    B.SetCartesianCoordinates(1., 2., 3.);
    EXPECT_TRUE(A != B);
    Vector3D C;
    C.SetSphericalCoordinates(0.1, 0.2, 0.3);
    EXPECT_TRUE(A != C);
}

TEST(Assignment, Copyconstructor)
{
    Vector3D A(1, 2, 3);
    Vector3D B(A);
    EXPECT_TRUE(A == B);
    Vector3D C(2, 3, 4);
    B = C;
    EXPECT_TRUE(C == B);
}

TEST(Assignment, Operator)
{
    Vector3D A;
    Vector3D B;
    A.SetCartesianCoordinates(1, 2, 3);
    EXPECT_TRUE(A != B);
    B = A;
    EXPECT_TRUE(A == B);
}

TEST(Assignment, Swap)
{
    Vector3D A;
    Vector3D B;
    EXPECT_TRUE(A == B);
    Vector3D* C = new Vector3D(1, 2, 3);
    Vector3D* D = new Vector3D(1, 2, 3);
    C->SetSphericalCoordinates(0.1, 0.2, 0.3);
    D->SetSphericalCoordinates(0.1, 0.2, 0.3);
    EXPECT_TRUE(*C == *D);

    A.swap(*C);
    EXPECT_TRUE(A == *D);
    EXPECT_TRUE(B == *C);
}

TEST(Addition, Operator)
{
    Vector3D A;
    Vector3D B;
    Vector3D C;
    Vector3D D;
    A.SetCartesianCoordinates(1, 2, 3);
    B.SetCartesianCoordinates(2, 4, 6);
    C.SetCartesianCoordinates(3, 6, 9);
    D = A + B;
    EXPECT_TRUE(C == D);
    D.SetCartesianCoordinates(0, 0, 0);
    EXPECT_TRUE(D != C);
    D = B + A;
    EXPECT_TRUE(D == C);
}

TEST(Subtraction, Operator)
{
    Vector3D A;
    Vector3D B;
    Vector3D C;
    Vector3D D;
    Vector3D E;
    A.SetCartesianCoordinates(1, 2, 3);
    B.SetCartesianCoordinates(2, 4, 6);
    C.SetCartesianCoordinates(3, 6, 9);
    E.SetCartesianCoordinates(-1, -2, -3);
    D = -A;
    EXPECT_TRUE(D == E);
    D = C - B;
    EXPECT_TRUE(D == A);
    D = B - C;
    EXPECT_TRUE(D == E);
}

TEST(Scaling, Operator)
{
    Vector3D A;
    Vector3D B;
    Vector3D C;
    double factor1 = 2.;
    A.SetCartesianCoordinates(1, 2, 3);
    B.SetCartesianCoordinates(2, 4, 6);
    C = A * factor1;
    EXPECT_TRUE(C == B);
    C.SetCartesianCoordinates(0, 0, 0);
    EXPECT_TRUE(C != B);
    C = factor1 * A;
    EXPECT_TRUE(C == B);
}

TEST(ScalarProduct, Operator)
{
    Vector3D A;
    Vector3D B;
    Vector3D C;
    double factor1 = 5;
    double result  = 0;
    A.SetCartesianCoordinates(1, 1, 1);
    B.SetCartesianCoordinates(2, 1, 2);
    EXPECT_TRUE(result != factor1);
    result = scalar_product(A, B);
    EXPECT_TRUE(result == factor1);
    result = 0.;
    result = scalar_product(B, A);
    EXPECT_TRUE(result == factor1);
}

TEST(VectorProduct, Operator)
{
    Vector3D A;
    Vector3D B;
    Vector3D C;
    Vector3D D;
    Vector3D E;
    A.SetCartesianCoordinates(1, 0, 0);
    B.SetCartesianCoordinates(0, 1, 0);
    C.SetCartesianCoordinates(0, 0, 1);
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
    Vector3D A;
    Vector3D B;
    double sum = 3;
    double result;
    A.SetCartesianCoordinates(1, 2, 2);
    B.SetCartesianCoordinates(1, 1, 1);
    result = A.magnitude();
    EXPECT_TRUE(result == sum);
    result = B.magnitude();
    EXPECT_TRUE(result != sum);
}

TEST(Normalise, Operator)
{
    Vector3D A;
    Vector3D B;
    Vector3D C;
    A.SetCartesianCoordinates(3, 0, 0);
    B.SetCartesianCoordinates(1, 1, 1);
    C.SetCartesianCoordinates(1, 0, 0);
    EXPECT_TRUE(A != C);
    EXPECT_TRUE(B != C);
    A.normalise();
    EXPECT_TRUE(A == C);
    B.normalise();
    EXPECT_TRUE(B != C);
}

TEST(CalculateSphericalCoordinates, Conversion)
{
    Vector3D A;
    Vector3D B(1, 2, 2);
    A.SetCartesianCoordinates(1, 2, 2);
    EXPECT_TRUE(A == B);
    A.CalculateSphericalCoordinates();
    EXPECT_TRUE(A != B);
    B.SetSphericalCoordinates(3, std::atan2(2., 1.), std::acos(2. / 3.));
    EXPECT_TRUE(B == A);
}

TEST(CalculateCartesianFromSpherical, Conversion)
{
    Vector3D A;
    Vector3D B;
    double epsilon = std::numeric_limits<double>::epsilon();
    A.SetCartesianCoordinates(1, 2, 2);
    B.SetSphericalCoordinates(3, std::atan2(2., 1.), std::acos(2. / 3.));
    B.CalculateCartesianFromSpherical();
    EXPECT_TRUE(A != B);
    B.SetSphericalCoordinates(0, 0, 0);
    double error_factor = 2.;
    bool test_x = std::abs(A.GetX() - B.GetX()) < std::abs(std::min(A.GetX(), B.GetX())) * epsilon * error_factor;
    bool test_y = std::abs(A.GetY() - B.GetY()) < std::abs(std::min(A.GetY(), B.GetY())) * epsilon * error_factor;
    bool test_z = std::abs(A.GetZ() - B.GetZ()) < std::abs(std::min(A.GetZ(), B.GetZ())) * epsilon * error_factor;
    EXPECT_TRUE(A == B || test_x && test_y && test_z);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
