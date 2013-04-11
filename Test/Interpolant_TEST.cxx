#include "gtest/gtest.h"
#include "PROPOSAL/Interpolant.h"
#include "math.h"

double X2(double x){
    return x*x*exp(x);
}

double X_YY(double x, double y){
    return x + y*y*exp(x);
}

int max = 100;
double xmin = 3;
double xmax = 20;
int romberg = 5;
bool rational = false;
bool relative = false;
bool isLog = false;
int rombergY = 5;
bool rationalY = false;
bool relativeY = false;
bool logSubst = false;
int max2 = 100;
double x2min = 5;
double x2max = 20;
int romberg2 = 5;
bool rational2 = false;
bool relative2 = false;
bool isLog2 = false;

TEST(_1D_Interpol , Simple_Test_of_XX_EXPX ) {
    Interpolant* Pol1 = new Interpolant(max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);
    double SearchX = 5;
    double precision = 1E-5;

    double PolValue = Pol1->interpolate(SearchX);
    double RealValue = X2(SearchX);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    delete Pol1;
}


TEST(_1D_Interpol , Save_To_File ) {
    Interpolant* Pol1 = new Interpolant(max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);

    ASSERT_TRUE(Pol1->Save(".TestSaveInterpol.txt"));

    delete Pol1;
}

TEST(_1D_Interpol , Load_From_File ) {
    Interpolant* Pol1 = new Interpolant();
    Pol1->Load(".TestSaveInterpol.txt");
    double SearchX = 5;
    double precision = 1E-5;

    double PolValue = Pol1->interpolate(SearchX);
    double RealValue = X2(SearchX);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    delete Pol1;
}

TEST(_1D_Interpol , Rational_On ) {
    Interpolant* Pol1 = new Interpolant(max, xmin, xmax, X2, romberg, true, relative, isLog, rombergY, rationalY, relativeY, logSubst);
    double SearchX = 5;
    double precision = 1E-5;

    double PolValue = Pol1->interpolate(SearchX);
    double RealValue = X2(SearchX);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    PolValue = Pol1->findLimit(RealValue);
    ASSERT_NEAR(PolValue, SearchX, SearchX*precision );

    delete Pol1;
}

TEST(_1D_Interpol , isLog_On ) {
    Interpolant* Pol1 = new Interpolant(max, xmin, xmax, X2, romberg, rational, relative, true, rombergY, rationalY, relativeY, logSubst);
    double SearchX = 5;
    double precision = 1E-5;

    double PolValue = Pol1->interpolate(SearchX);
    double RealValue = X2(SearchX);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    PolValue = Pol1->findLimit(RealValue);
    ASSERT_NEAR(PolValue, SearchX, SearchX*precision );

    delete Pol1;
}

TEST(_1D_Interpol , RationalY_On ) {
    Interpolant* Pol1 = new Interpolant(max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, !rationalY, relativeY, logSubst);
    double SearchX = 5;
    double precision = 1E-4;

    double PolValue = Pol1->interpolate(SearchX);
    double RealValue = X2(SearchX);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    PolValue = Pol1->findLimit(RealValue);
    ASSERT_NEAR(PolValue, SearchX, SearchX*precision );

    delete Pol1;
}

TEST(_1D_Interpol , logSubst_On ) {
    Interpolant* Pol1 = new Interpolant(max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, !logSubst);
    double SearchX = 5;
    double precision = 1E-5;

    double PolValue = Pol1->interpolate(SearchX);
    double RealValue = X2(SearchX);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    PolValue = Pol1->findLimit(RealValue);
    ASSERT_NEAR(PolValue, SearchX, SearchX*precision );

    delete Pol1;
}

TEST(_1D_Interpol , All_On ) {
    Interpolant* Pol1 = new Interpolant(max, xmin, xmax, X2, romberg, !rational, relative, !isLog, rombergY, !rationalY, relativeY, !logSubst);
    double SearchX = 5;
    double precision = 1E-5;

    double PolValue = Pol1->interpolate(SearchX);
    double RealValue = X2(SearchX);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    PolValue = Pol1->findLimit(RealValue);
    ASSERT_NEAR(PolValue, SearchX, SearchX*precision );

    delete Pol1;
}


TEST(_2D_Interpol , Simple_Test_of_X_YY_EXPX ) {
    Interpolant* Pol2 = new Interpolant(max,  xmin, xmax, max2, x2min, x2max, X_YY,
                                        romberg, rational, relative, isLog,
                                        romberg2, rational2, relative2, isLog2,
                                        rombergY, rationalY, relativeY, logSubst);
    double SearchX = 7;
    double SearchY = 11;
    double precision = 1E-6;

    double PolValue = Pol2->interpolate(SearchX,SearchY);
    double RealValue = X_YY(SearchX,SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    delete Pol2;
}

TEST(_2D_Interpol , rational1_On ) {
    Interpolant* Pol2 = new Interpolant(max,  xmin, xmax, max2, x2min, x2max, X_YY,
                                        romberg, !rational, relative, isLog,
                                        romberg2, rational2, relative2, isLog2,
                                        rombergY, rationalY, relativeY, logSubst);
    double SearchX = 7;
    double SearchY = 11;
    double precision = 1E-6;

    double PolValue = Pol2->interpolate(SearchX,SearchY);
    double RealValue = X_YY(SearchX,SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    PolValue = Pol2->findLimit(SearchX,RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY*precision );

    delete Pol2;
}

TEST(_2D_Interpol , rational2_On ) {
    Interpolant* Pol2 = new Interpolant(max,  xmin, xmax, max2, x2min, x2max, X_YY,
                                        romberg, rational, relative, isLog,
                                        romberg2, !rational2, relative2, isLog2,
                                        rombergY, rationalY, relativeY, logSubst);
    double SearchX = 7;
    double SearchY = 11;
    double precision = 1E-6;

    double PolValue = Pol2->interpolate(SearchX,SearchY);
    double RealValue = X_YY(SearchX,SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    PolValue = Pol2->findLimit(SearchX,RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY*precision );

    delete Pol2;
}

TEST(_2D_Interpol , rational12_On ) {
    Interpolant* Pol2 = new Interpolant(max,  xmin, xmax, max2, x2min, x2max, X_YY,
                                        romberg, !rational, relative, isLog,
                                        romberg2, !rational2, relative2, isLog2,
                                        rombergY, rationalY, relativeY, logSubst);
    double SearchX = 7;
    double SearchY = 11;
    double precision = 1E-6;

    double PolValue = Pol2->interpolate(SearchX,SearchY);
    double RealValue = X_YY(SearchX,SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    PolValue = Pol2->findLimit(SearchX,RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY*precision );

    delete Pol2;
}

TEST(_2D_Interpol , isLog1_On ) {
    Interpolant* Pol2 = new Interpolant(max,  xmin, xmax, max2, x2min, x2max, X_YY,
                                        romberg, rational, relative, !isLog,
                                        romberg2, rational2, relative2, isLog2,
                                        rombergY, rationalY, relativeY, logSubst);
    double SearchX = 7;
    double SearchY = 11;
    double precision = 1E-6;

    double PolValue = Pol2->interpolate(SearchX,SearchY);
    double RealValue = X_YY(SearchX,SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    PolValue = Pol2->findLimit(SearchX,RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY*precision );

    delete Pol2;
}

TEST(_2D_Interpol , isLog2_On ) {
    Interpolant* Pol2 = new Interpolant(max,  xmin, xmax, max2, x2min, x2max, X_YY,
                                        romberg, rational, relative, isLog,
                                        romberg2, rational2, relative2, !isLog2,
                                        rombergY, rationalY, relativeY, logSubst);
    double SearchX = 7;
    double SearchY = 11;
    double precision = 1E-6;

    double PolValue = Pol2->interpolate(SearchX,SearchY);
    double RealValue = X_YY(SearchX,SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    PolValue = Pol2->findLimit(SearchX,RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY*precision );

    delete Pol2;
}

TEST(_2D_Interpol , isLog12_On ) {
    Interpolant* Pol2 = new Interpolant(max,  xmin, xmax, max2, x2min, x2max, X_YY,
                                        romberg, rational, relative, !isLog,
                                        romberg2, rational2, relative2, !isLog2,
                                        rombergY, rationalY, relativeY, logSubst);
    double SearchX = 7;
    double SearchY = 11;
    double precision = 1E-6;

    double PolValue = Pol2->interpolate(SearchX,SearchY);
    double RealValue = X_YY(SearchX,SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    PolValue = Pol2->findLimit(SearchX,RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY*precision );

    delete Pol2;
}

TEST(_2D_Interpol , rationalY_On ) {
    Interpolant* Pol2 = new Interpolant(max,  xmin, xmax, max2, x2min, x2max, X_YY,
                                        romberg, rational, relative, isLog,
                                        romberg2, rational2, relative2, isLog2,
                                        rombergY, !rationalY, relativeY, logSubst);
    double SearchX = 7;
    double SearchY = 11;
    double precision = 1E-6;

    double PolValue = Pol2->interpolate(SearchX,SearchY);
    double RealValue = X_YY(SearchX,SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    PolValue = Pol2->findLimit(SearchX,RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY*precision );

    delete Pol2;
}

TEST(_2D_Interpol , logSubst_On ) {
    Interpolant* Pol2 = new Interpolant(max,  xmin, xmax, max2, x2min, x2max, X_YY,
                                        romberg, rational, relative, isLog,
                                        romberg2, rational2, relative2, isLog2,
                                        rombergY, rationalY, relativeY, !logSubst);
    double SearchX = 7;
    double SearchY = 11;
    double precision = 1E-6;

    double PolValue = Pol2->interpolate(SearchX,SearchY);
    double RealValue = X_YY(SearchX,SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    PolValue = Pol2->findLimit(SearchX,RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY*precision );

    delete Pol2;
}

TEST(_2D_Interpol , All_On ) {
    Interpolant* Pol2 = new Interpolant(max,  xmin, xmax, max2, x2min, x2max, X_YY,
                                        romberg, !rational, relative, !isLog,
                                        romberg2, !rational2, relative2, !isLog2,
                                        rombergY, !rationalY, relativeY, !logSubst);
    double SearchX = 7;
    double SearchY = 11;
    double precision = 1E-6;

    double PolValue = Pol2->interpolate(SearchX,SearchY);
    double RealValue = X_YY(SearchX,SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    PolValue = Pol2->findLimit(SearchX,RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY*precision );

    delete Pol2;
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
