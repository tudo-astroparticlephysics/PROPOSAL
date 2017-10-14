
// #include <iostream>
// #include <string>
#include <cmath>

#include "gtest/gtest.h"

#include "PROPOSAL/math/RootFinder.h"

using namespace PROPOSAL;

double Exp(double r){
    return exp(r) - exp(0);
}

double DiffExp(double r){
    return exp(r);
}

double Polynom(double r){
    return r*r;
}

double DiffPolynom(double r){
    return 2*r - 4;
}

double Polynom2(double r){
    return r*r -4;
}

double DiffPolynom2(double r){
    return 2*r;
}

double Polynom3(double r){
    return -r*r -4;
}

double DiffPolynom3(double r){
    return -2*r;
}


TEST(Comparison , Comparison_equal ) {
    RootFinder A;
    RootFinder B;
    EXPECT_TRUE(A==B);
    RootFinder* C = new RootFinder(10,1e-3);
    RootFinder* D = new RootFinder(10,1e-3);
    EXPECT_TRUE(*C==*D);
    RootFinder* E = new RootFinder(20,1e-6);
    EXPECT_TRUE(A==*E);

}

TEST(Comparison , Comparison_not_equal ) {
    RootFinder A;
    RootFinder B(10,1e-3);
    EXPECT_TRUE(A!=B);
    RootFinder* C = new RootFinder(10,1e-3);
    RootFinder* D = new RootFinder(20,0.3);
    EXPECT_TRUE(*C!=*D);
    RootFinder* E = new RootFinder(10,1e-3);
    EXPECT_TRUE(*C==*E);
    E->SetMaxSteps(20);
    EXPECT_TRUE(*C!=*E);


}

TEST(Assignment , Copyconstructor ) {
    RootFinder A;
    RootFinder B =A;

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Copyconstructor2 ) {
    RootFinder A(10,1e-3);
    RootFinder B(A);

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Operator ) {
    RootFinder A;
    RootFinder B(10,0.4);

    EXPECT_TRUE(A!=B);

    B=A;

    EXPECT_TRUE(A==B);
}

TEST(Assignment , Swap ) {
    RootFinder A;
    RootFinder B;
    EXPECT_TRUE(A==B);
    RootFinder* C = new RootFinder(5,0.3);
    RootFinder* D = new RootFinder(5,0.3);
    EXPECT_TRUE(*C==*D);

    A.swap(*C);
    EXPECT_TRUE(A==*D);
    EXPECT_TRUE(B==*C);


}


TEST(RootFinder , e_to_x_minus_e ) {
    RootFinder *finder = new RootFinder();
    ASSERT_NEAR(finder->FindRoot(-2,2,1, Exp, DiffExp) ,0, 1E-10);
    delete finder;
}

TEST(RootFinder , x_times_x_intersection_with_4 ) {
    RootFinder *finder = new RootFinder();
    ASSERT_NEAR(finder->FindRoot(-2,2,1, Polynom, DiffPolynom) ,-2., 1E-6*0);
    delete finder;
}

TEST(RootFinder , x_times_x_minus_4 ) {
    RootFinder *finder = new RootFinder();
    ASSERT_NEAR(finder->FindRoot(-3,3,1, Polynom2, DiffPolynom2) ,-2, 0.);
    delete finder;
}

TEST(RootFinder , minus_x_times_x_intersection_with_minus4 ) {
    RootFinder *finder = new RootFinder();
    ASSERT_NEAR(finder->FindRoot(-2,2,0, Polynom3, DiffPolynom3) ,-2., 1E-6*0);
    delete finder;
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

//TEST(IntegralValue, EqualBorders) {
//    Integral* Int = new Integral();

//    EXPECT_EQ(Int->Integrate(3,3,Testfkt,1),0);

//    delete Int;
//}

//TEST(IntegralValue, SmallError) {
//    Integral* Int = new Integral();

//    ASSERT_NEAR(   Int->Integrate(0,3,Testexp,1),exp(3)-1
//                            ,(exp(3)-1)*1.e-6);

//    delete Int;
//}

//TEST(IntegralValue, FloatEqual) {
//    Integral* Int = new Integral();
//    //Last 4 digits can differ. relError<1E-4
//    ASSERT_FLOAT_EQ(Int->Integrate(0,3,Testexp,1),exp(3)-1);

//    delete Int;
//}
