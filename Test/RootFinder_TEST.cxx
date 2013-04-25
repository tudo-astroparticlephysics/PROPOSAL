#include "gtest/gtest.h"
#include "PROPOSAL/RootFinder.h"
#include <iostream>
#include <string>
#include <cmath>

double Exp(double r){
    return exp(r) - exp(0);
}

double DiffExp(double r){
    return exp(r);
}

TEST(RootFinder , e_to_x_minus_e ) {
    RootFinder *finder = new RootFinder();
    ASSERT_NEAR(finder->FindRoot(-2,2,1, Exp, DiffExp,0.0001) ,0, 1E-6);
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
