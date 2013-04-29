#include "gtest/gtest.h"
#include <limits.h>
#include "PROPOSAL/Integral.h"
#include "math.h"
#include <iostream>

double Testfkt(double r){
  return exp(r);
}

double Testexp(double r){
    return exp(r);
}

bool relErr(double Is, double HasToBe, double RelError){
    return fabs((Is - HasToBe)/HasToBe)<RelError;
}

TEST(Assignment , Copyconstructor ) {
    Integral A;
    Integral B =A;

    EXPECT_EQ(A.GetMaxSteps() , B.GetMaxSteps());
    EXPECT_EQ(A.GetRomberg()  , B.GetRomberg());
    EXPECT_EQ(A.GetIX().size(), B.GetIX().size());
    EXPECT_EQ(A.GetIY().size(), B.GetIY().size());
    EXPECT_EQ(A.GetC().size() , B.GetC().size());
    EXPECT_EQ(A.GetD().size() , B.GetD().size());

    for(int i =0; i < A.GetIX().size(); i++){
        EXPECT_EQ(A.GetIX().at(i), B.GetIX().at(i));
    }
    for(int i =0; i < A.GetIY().size(); i++){
        EXPECT_EQ(A.GetIY().at(i), B.GetIY().at(i));
    }
    for(int i =0; i < A.GetC().size(); i++){
        EXPECT_EQ(A.GetC().at(i), B.GetC().at(i));
    }
    for(int i =0; i < A.GetD().size(); i++){
        EXPECT_EQ(A.GetD().at(i), B.GetD().at(i));
    }

}

TEST(Assignment , Copyconstructor2 ) {
    Integral A;
    Integral B(A);

    EXPECT_EQ(A.GetMaxSteps() , B.GetMaxSteps());
    EXPECT_EQ(A.GetRomberg()  , B.GetRomberg());
    EXPECT_EQ(A.GetIX().size(), B.GetIX().size());
    EXPECT_EQ(A.GetIY().size(), B.GetIY().size());
    EXPECT_EQ(A.GetC().size() , B.GetC().size());
    EXPECT_EQ(A.GetD().size() , B.GetD().size());

    for(int i =0; i < A.GetIX().size(); i++){
        EXPECT_EQ(A.GetIX().at(i), B.GetIX().at(i));
    }
    for(int i =0; i < A.GetIY().size(); i++){
        EXPECT_EQ(A.GetIY().at(i), B.GetIY().at(i));
    }
    for(int i =0; i < A.GetC().size(); i++){
        EXPECT_EQ(A.GetC().at(i), B.GetC().at(i));
    }
    for(int i =0; i < A.GetD().size(); i++){
        EXPECT_EQ(A.GetD().at(i), B.GetD().at(i));
    }
}



TEST(Assignment , Operator ) {
    Integral A;
    A.Integrate(0,3,Testfkt,1);
    Integral B(8,40,1e-9);

    EXPECT_NE(A.GetMaxSteps() , B.GetMaxSteps());
    EXPECT_NE(A.GetRomberg()  , B.GetRomberg());
    EXPECT_NE(A.GetIX().size(), B.GetIX().size());
    EXPECT_NE(A.GetIY().size(), B.GetIY().size());
    EXPECT_NE(A.GetC().size() , B.GetC().size());
    EXPECT_NE(A.GetD().size() , B.GetD().size());

    B=A;

    EXPECT_EQ(A.GetMaxSteps() , B.GetMaxSteps());
    EXPECT_EQ(A.GetRomberg()  , B.GetRomberg());
    EXPECT_EQ(A.GetIX().size(), B.GetIX().size());
    EXPECT_EQ(A.GetIY().size(), B.GetIY().size());
    EXPECT_EQ(A.GetC().size() , B.GetC().size());
    EXPECT_EQ(A.GetD().size() , B.GetD().size());

    for(int i =0; i < A.GetIX().size(); i++){
        EXPECT_EQ(A.GetIX().at(i), B.GetIX().at(i));
    }
    for(int i =0; i < A.GetIY().size(); i++){
        EXPECT_EQ(A.GetIY().at(i), B.GetIY().at(i));
    }
    for(int i =0; i < A.GetC().size(); i++){
        EXPECT_EQ(A.GetC().at(i), B.GetC().at(i));
    }
    for(int i =0; i < A.GetD().size(); i++){
        EXPECT_EQ(A.GetD().at(i), B.GetD().at(i));
    }
}

TEST(IntegralValue , Zero_to_Three_of_xx ) {
    Integral* Int = new Integral();
    ASSERT_NEAR(Int->Integrate(0,3,Testfkt,1),exp(3)-1 , (exp(3)-1)*1E-6);
    delete Int;
}

TEST(IntegralValue, EqualBorders) {
    Integral* Int = new Integral();

    EXPECT_EQ(Int->Integrate(3,3,Testfkt,1),0);

    delete Int;
}

TEST(IntegralValue, SmallError) {
    Integral* Int = new Integral();

    ASSERT_NEAR(   Int->Integrate(0,3,Testexp,1),exp(3)-1
                            ,(exp(3)-1)*1.e-6);

    delete Int;
}

TEST(IntegralValue, FloatEqual) {
    Integral* Int = new Integral();
    //Last 4 digits can differ. relError<1E-4
    ASSERT_FLOAT_EQ(Int->Integrate(0,3,Testexp,1),exp(3)-1);

    delete Int;
}

TEST(IntegralValue, MultiplePrecisions) {
    double xmin=0,xmax=3;
    double  ExactIntegral=exp(3)-1;
    double  CalcIntegral=0;

    double precision = 1E-5;
    for(double precision = 1E-5; precision>1E-16;precision/=10){

        Integral* Int = new Integral(5,20,precision);
        CalcIntegral = Int->Integrate(xmin,xmax,Testexp,1);

        ASSERT_NEAR(CalcIntegral,ExactIntegral, ExactIntegral*precision);

        delete Int;
    }
}

TEST(IntegralValue, IntegrateWithSubstitution) {
    double xmin=2,xmax=4;
    double  ExactIntegral=exp(xmax)-exp(xmin);
    double  CalcIntegral=0;

    double precision = 1E-5;
    for(double precision = 1E-5; precision>1E-11;precision/=10){

        Integral* Int = new Integral(5,20,precision);
        CalcIntegral = Int->Integrate(xmin,xmax,Testexp,3,2.);

        ASSERT_NEAR(CalcIntegral,ExactIntegral, ExactIntegral*precision);

        delete Int;
    }
}

TEST(IntegralValue, IntegrateWithLog) {
    double xmin=2,xmax=4;
    double  ExactIntegral=exp(xmax)-exp(xmin);
    double  CalcIntegral=0;

    double precision = 1E-5;
    for(double precision = 1E-5; precision>1E-11;precision/=10){

        Integral* Int = new Integral(5,20,precision);
        CalcIntegral = Int->Integrate(xmin,xmax,Testexp,4);

        ASSERT_NEAR(CalcIntegral,ExactIntegral, ExactIntegral*precision);

        delete Int;
    }
}

TEST(IntegralValue, IntegrateWithLogSubstitution) {
    double xmin=2,xmax=4;
    double  ExactIntegral=exp(xmax)-exp(xmin);
    double  CalcIntegral=0;

    double precision = 1E-5;
    for(double precision = 1E-5; precision>1E-6;precision/=10){

        Integral* Int = new Integral(5,20,precision);
        CalcIntegral = Int->Integrate(xmin,xmax,Testexp,5,2.);

        ASSERT_NEAR(CalcIntegral,ExactIntegral, ExactIntegral*precision);

        delete Int;
    }
}

//TEST(IntegralValue, HasToFail) {
//    EXPECT_TRUE(false);
//}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
