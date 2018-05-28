
#include "gtest/gtest.h"

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"

using namespace PROPOSAL;

TEST(Comparison, Comparison_equal)
{
    Medium* A = new AntaresWater(0.3);
    Medium* B = new AntaresWater(0.3);
    EXPECT_TRUE(*A == *B);

    Components::Oxygen a;
    Components::Oxygen b;
    EXPECT_TRUE(a == b);

    Medium* C = new Water();
    Water D;
    EXPECT_TRUE(D == *C);

    Medium* E = MediumFactory::Get().CreateMedium("IcE");
    Ice F;
    EXPECT_TRUE(F == *E);

    Water G;
    Water H;
    EXPECT_TRUE(G == H);

    delete A;
    delete B;
    delete C;
    delete E;
}

TEST(Comparison, Comparison_not_equal)
{
    AntaresWater A;
    AntaresWater B(0.3);
    EXPECT_TRUE(A != B);

    Medium* C = new AntaresWater(0.3);
    Medium* D = new Ice(0.3);
    EXPECT_TRUE(*C != *D);

    Medium* E = new Water(0.3);
    Medium* F = new Water(1.0);
    EXPECT_TRUE(*E != *F);

    Medium* G = MediumFactory::Get().CreateMedium("WaTeR");
    EXPECT_TRUE(*E != *G);

    Components::Hydrogen a;
    Components::Oxygen b;
    EXPECT_TRUE(a != b);

    Components::Component* c = new Components::Oxygen();
    Components::Component* d = new Components::Oxygen(2.0);
    EXPECT_TRUE(c != d);

    Components::Component* e = new Components::Magnesium();
    Components::Component* f = new Components::Iron();
    EXPECT_TRUE(e != f);

    delete c;
    delete d;
    delete e;
    delete f;

    delete C;
    delete D;
    delete E;
    delete F;
    delete G;
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

    Components::Component* c = new Components::Hydrogen();
    Components::Component* d = new Components::Hydrogen(2.0);
    *c                       = *d;

    EXPECT_TRUE(*c == *d);

    Water A;
    Water B;
    B = A;

    EXPECT_TRUE(A == B);

    Medium* C = new Water();
    Medium* D = new Water(2.0);
    *D        = *C;

    EXPECT_TRUE(*D == *C);

    Medium* E = new Water();
    Medium* F = new Ice();
    *E        = *F;

    EXPECT_TRUE(*E != *F);

    delete c;
    delete d;
    delete C;
    delete D;
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

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
