#include "gtest/gtest.h"
#include "PROPOSAL/Integral.h"

double Testfkt(double r){
  return r*r;
}


TEST(FooTest, DoesXyz) {
    Integral* Int = new Integral();
    EXPECT_EQ(Int->IntegrateClosed(0,3,Testfkt),9);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}