
#include "gtest/gtest.h"

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/medium/density_distr/density_exponential.h"
#include "PROPOSAL/medium/density_distr/density_homogeneous.h"
#include "PROPOSAL/medium/density_distr/density_polynomial.h"
#include "PROPOSAL/medium/density_distr/density_splines.h"

using namespace PROPOSAL;

TEST(Comparison, Comparison_equal)
{
    Axis* ax_A = new RadialAxis();
    Axis* ax_B = new RadialAxis();
    EXPECT_TRUE(*ax_A == *ax_B);

    Axis* ax_C = new CartesianAxis();
    CartesianAxis ax_D;
    EXPECT_TRUE(*ax_C == ax_D);

    CartesianAxis ax_E;
    EXPECT_TRUE(ax_E == ax_D);

    Density_distr* A = new Density_homogeneous(0.3);
    Density_distr* B = new Density_homogeneous(0.3);
    EXPECT_TRUE(*A == *B);

    double sigma = 1.;
    Density_distr* C = new Density_exponential(ax_D, sigma);
    Density_exponential D(ax_D, sigma);
    EXPECT_TRUE(D == *C);

    Density_homogeneous G;
    Density_homogeneous H;
    EXPECT_TRUE(G == H);

    delete A;
    delete B;
    delete C;
    delete ax_A;
    delete ax_B;
    delete ax_C;
}

TEST(Comparison, Comparison_not_equal)
{
    Vector3D faxis_default(1,0,0);
    Vector3D faxis_new(0,1,0);
    Vector3D fp0_default(0,0,0);
    Vector3D fp0_new(1,1,1);
    CartesianAxis ax_A;
    CartesianAxis ax_B(faxis_default, fp0_new);
    CartesianAxis ax_C(faxis_new, fp0_default);
    EXPECT_TRUE(ax_A != ax_B);
    EXPECT_TRUE(ax_A != ax_C);

    Axis* ax_D = new CartesianAxis();
    Axis* ax_E = new RadialAxis();
    EXPECT_TRUE(*ax_D != *ax_E);

    Density_distr* A = new Density_homogeneous();
    Density_distr* B = new Density_exponential(ax_A, 1);
    Density_distr* D = new Density_exponential(ax_A, 2);
    EXPECT_TRUE(*A != *B);
    EXPECT_TRUE(*B != *D);

    std::vector<double> vecA = {1,2};
    std::vector<double> vecB = {2,3};
    Polynom poly_A(vecA);
    Polynom poly_B(vecB);
    Density_polynomial E(ax_A, poly_A);
    Density_polynomial F(ax_B, poly_A);
    Density_polynomial G(ax_A, poly_B);
    EXPECT_TRUE(E != F);
    EXPECT_TRUE(E != G);

    Linear_Spline spl_A(vecA, vecA);
    Linear_Spline spl_B(vecA, vecB);
    Density_splines H(ax_A, spl_A);
    Density_splines I(ax_B, spl_A);
    Density_splines J(ax_A, spl_B);
    EXPECT_TRUE(H != I);
    EXPECT_TRUE(H != J);

    delete ax_D;
    delete ax_E;
    delete A;
    delete B;
    delete D;
}

TEST(Assignment, Copyconstructor)
{
    CartesianAxis ax_A;
    CartesianAxis ax_B = ax_A;
    CartesianAxis ax_C(ax_A);
    EXPECT_TRUE(ax_A == ax_B);
    EXPECT_TRUE(ax_A == ax_C);

    std::vector<double> vecA = {1,2};
    Polynom poly_A(vecA);
    Density_polynomial A(ax_A, poly_A);
    Density_polynomial B = A;
    Density_polynomial C(A);
    EXPECT_TRUE(A == B);
    EXPECT_TRUE(A == C);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
