#include "gtest/gtest.h"
#include "PROPOSAL/Interpolant.h"
#include "math.h"

double X2(double x){
    return x*x;
}

double X_YY(double x, double y){
    return x + y*y;
}

int max = 100;
double xmin = 0;
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
double x2min = 0;
double x2max = 20;
int romberg2 = 5;
bool rational2 = false;
bool relative2 = false;
bool isLog2 = false;

TEST(_1D_Interpol , Test_of_XX ) {
    Interpolant* Pol1 = new Interpolant(max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);
    double SearchX = 3;
    double precision = 1E-6;

    double PolValue = Pol1->interpolate(SearchX);
    double RealValue = X2(SearchX);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    delete Pol1;
}


TEST(_2D_Interpol , Test_of_X_YY ) {
    Interpolant* Pol2 = new Interpolant(max,  xmin, xmax, max2, x2min, x2max, X_YY,
                                        romberg, rational, relative, isLog,
                                        romberg2, rational2, relative2, isLog2,
                                        rombergY, rationalY, relativeY, logSubst);
    double SearchX = 2.;
    double SearchY = 3.;
    double precision = 1E-6;

    double PolValue = Pol2->interpolate(SearchX,SearchY);
    double RealValue = X_YY(SearchX,SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue*precision );

    delete Pol2;
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
