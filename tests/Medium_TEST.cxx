
#include "gtest/gtest.h"

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"

#include "PROPOSAL/density_distr/density_homogeneous.h"

using namespace PROPOSAL;

TEST(Comparison, Comparison_equal)
{
    std::unique_ptr<Medium> A(new AntaresWater());
    std::unique_ptr<Medium> B(new AntaresWater());
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
    A.SetMassDensity(2. * A.GetMassDensity());
    AntaresWater B;
    EXPECT_TRUE(A != B);

    std::unique_ptr<Medium> C(new AntaresWater());
    std::unique_ptr<Medium> D(new Ice());
    EXPECT_TRUE(*C != *D);

    std::unique_ptr<Medium> E(new Water());
    E->SetMassDensity(2. * E->GetMassDensity());
    std::unique_ptr<Medium> F(new Water());
    EXPECT_TRUE(*E != *F);

    std::shared_ptr<const Medium> G = CreateMedium("WaTeR");
    EXPECT_TRUE(*E != *G);

    Components::Hydrogen a;
    Components::Oxygen b;
    EXPECT_TRUE(a != b);
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

    std::unique_ptr<Component> c(new Components::Hydrogen());
    std::unique_ptr<Component> d(new Components::Hydrogen());
    *c                       = *d;

    EXPECT_TRUE(*c == *d);

    Water A;
    Water B;
    B = A;

    EXPECT_TRUE(A == B);

    std::unique_ptr<Medium> C(new Water());
    std::unique_ptr<Medium> D(new Water());
    *D        = *C;

    EXPECT_TRUE(*D == *C);

    std::unique_ptr<Medium> E(new Water());
    std::unique_ptr<Medium> F(new Ice());
    *E        = *F;

    EXPECT_TRUE(*E == *F);
}

TEST(Assignment, Copyconstructor2)
{
    Uranium A;
    Uranium B(A);

    EXPECT_TRUE(A == B);
}

// TEST(Assignment, Swap)
// {
//     Water A;
//     Water B;
//     EXPECT_TRUE(A == B);
//     Medium* C = new AntaresWater();
//     Medium* D = new AntaresWater();
//     EXPECT_TRUE(*C == *D);
//     Medium* E = new Water();
//     EXPECT_TRUE(A == *E);

//     A.swap(*C);
//     EXPECT_TRUE(A == *D);
//     EXPECT_TRUE(B == *C);
//     A.swap(*C);
//     C->swap(*E);
//     EXPECT_TRUE(A == *C);

//     Medium* X = new Ice();
//     Medium* Y = new Copper();

//     X->swap(*Y);

//     Copper Z();
//     EXPECT_TRUE(*X == Z);
// }

// TEST(Density_distr, Evaluate)
// {
//     Water A;
//     Water B;

//     double corr_faktor = 1.0;
//     Density_homogeneous dens_hom(corr_faktor);
//     A.SetDensityDistribution(dens_hom);
//     B.SetDensityDistribution(dens_hom);

//     Vector3D test_point (2,1,0);

//     double DensA = A.GetDensityDistribution().Evaluate(test_point);
//     double DensB = B.GetDensityDistribution().Evaluate(test_point);

//     EXPECT_TRUE(DensA == DensB);
// }

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
