#include "gtest/gtest.h"
#include "PROPOSAL/Integral.h"
#include "math.h"

double Testfkt(double r){
  return r*r;
}

double Testexp(double r){
    return exp(r);
}

bool relErr(double Is, double HasToBe, double RelError){
    return fabs((Is - HasToBe)/HasToBe)<RelError;
}

TEST(IntegralValue , Zero_to_Three_of_xx ) {
    Integral* Int = new Integral();
    EXPECT_EQ(Int->IntegrateClosed(0,3,Testfkt),9);
}

TEST(IntegralValue, EqualBorders) {
    Integral* Int = new Integral();
    EXPECT_EQ(Int->IntegrateClosed(3,3,Testfkt),0);
}

TEST(IntegralValue, SmallError) {
    Integral* Int = new Integral();
    EXPECT_TRUE(  relErr(   Int->IntegrateClosed(0,3,Testexp),exp(3)-1
                            ,1.e-6) //Later GetPrecision();
                );
}

TEST(IntegralValue, HasToFail) {
    EXPECT_TRUE(false);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
