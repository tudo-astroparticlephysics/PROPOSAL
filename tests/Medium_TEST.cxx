
#include "gtest/gtest.h"

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"

#include "PROPOSAL/medium/density_distr/density_homogeneous.h"

using namespace PROPOSAL;

TEST(Comparison, Comparison_equal)
{
    std::unique_ptr<Medium> A(new AntaresWater(0.3));
    std::unique_ptr<Medium> B(new AntaresWater(0.3));
    EXPECT_TRUE(*A == *B);

    Components::Oxygen a;
    Components::Oxygen b;
    EXPECT_TRUE(a == b);

    std::unique_ptr<Medium> C(new Water());
    Water D;
    EXPECT_TRUE(D == *C);

    std::shared_ptr<const Medium> E = CreateMedium("IcE");
    Ice F;
    EXPECT_TRUE(F == *E);

    Water G;
    Water H;
    EXPECT_TRUE(G == H);
}

TEST(Comparison, Comparison_not_equal)
{
    AntaresWater A;
    AntaresWater B(0.3);
    EXPECT_TRUE(A != B);

    std::unique_ptr<Medium> C(new AntaresWater(0.3));
    std::unique_ptr<Medium> D(new Ice(0.3));
    EXPECT_TRUE(*C != *D);

    std::unique_ptr<Medium> E(new Water(0.3));
    std::unique_ptr<Medium> F(new Water(1.0));
    EXPECT_TRUE(*E != *F);

    std::shared_ptr<const Medium> G = CreateMedium("WaTeR");
    EXPECT_TRUE(*E != *G);

    Components::Hydrogen a;
    Components::Oxygen b;
    EXPECT_TRUE(a != b);

    std::unique_ptr<Components::Component> c(new Components::Oxygen());
    std::unique_ptr<Components::Component> d(new Components::Oxygen(2.0));
    EXPECT_TRUE(*c != *d);

    std::unique_ptr<Components::Component> e(new Components::Magnesium());
    std::unique_ptr<Components::Component> f(new Components::Iron());
    EXPECT_TRUE(*e != *f);
}

TEST(Assignment, Copyconstructor)
{
    Components::Hydrogen a;
    Components::Hydrogen b = a;

    EXPECT_TRUE(a == b);

    Water A;
    Water B = A;

    EXPECT_TRUE(A == B);
}

TEST(Assignment, AssignmentOperator)
{
    Components::Hydrogen a;
    Components::Hydrogen b;
    b = a;

    EXPECT_TRUE(a == b);

    std::unique_ptr<Components::Component> c(new Components::Hydrogen());
    std::unique_ptr<Components::Component> d(new Components::Hydrogen(2.0));
    *c                       = *d;

    EXPECT_TRUE(*c == *d);

    Water A;
    Water B;
    B = A;

    EXPECT_TRUE(A == B);

    std::unique_ptr<Medium> C(new Water());
    std::unique_ptr<Medium> D(new Water(2.0));
    *D        = *C;

    EXPECT_TRUE(*D == *C);

    std::unique_ptr<Medium> E(new Water());
    std::unique_ptr<Medium> F(new Ice());
    *E        = *F;

    EXPECT_TRUE(*E != *F);
}

TEST(Assignment, Copyconstructor2)
{
    Uranium A(1.3);
    Uranium B(A);

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Swap)
{
    Water A;
    Water B;
    EXPECT_TRUE(A == B);
    Medium* C = new AntaresWater(0.3);
    Medium* D = new AntaresWater(0.3);
    EXPECT_TRUE(*C == *D);
    Medium* E = new Water();
    EXPECT_TRUE(A == *E);

    A.swap(*C);
    EXPECT_TRUE(A == *D);
    EXPECT_TRUE(B == *C);
    A.swap(*C);
    C->swap(*E);
    EXPECT_TRUE(A == *C);

    Medium* X = new Ice(0.3);
    Medium* Y = new Copper(0.3);

    X->swap(*Y);

    Copper Z(0.3);
    EXPECT_TRUE(*X == Z);
}

TEST(Density_distr, Evaluate)
{
    Water A;
    Water B;

    double corr_faktor = 1.0;
    Density_homogeneous dens_hom(corr_faktor);
    A.SetDensityDistribution(dens_hom);
    B.SetDensityDistribution(dens_hom);

    Vector3D test_point (2,1,0);

    double DensA = A.GetDensityDistribution().Evaluate(test_point);
    double DensB = B.GetDensityDistribution().Evaluate(test_point);

    EXPECT_TRUE(DensA == DensB);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
