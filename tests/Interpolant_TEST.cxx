
#include <cmath>
#include "gtest/gtest.h"
#include "PROPOSAL/math/Interpolant.h"

using namespace PROPOSAL;

double X2(double x)
{
    return x * x * std::exp(x);
}

double X_YY(double x, double y)
{
    return x + y * y * std::exp(x);
}

int max        = 100;
double xmin    = 3;
double xmax    = 20;
int romberg    = 5;
bool rational  = false;
bool relative  = false;
bool isLog     = false;
int rombergY   = 5;
bool rationalY = false;
bool relativeY = false;
bool logSubst  = false;
int max2       = 100;
double x2min   = 5;
double x2max   = 20;
int romberg2   = 5;
bool rational2 = false;
bool relative2 = false;
bool isLog2    = false;

std::string File1DTest = "Interpol1D_Save.txt";
std::string File2DTest = "Interpol2D_Save.txt";

TEST(Comparison, Comparison_equal)
{

    double SearchX = 5;
    double PolValue;
    Interpolant A;
    Interpolant B;
    EXPECT_TRUE(A == B);
    Interpolant* C = new Interpolant(
        max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);
    Interpolant* D = new Interpolant(
        max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);
    EXPECT_TRUE(*C == *D);
    PolValue = C->Interpolate(SearchX);
    PolValue = D->Interpolate(SearchX);
    EXPECT_TRUE(*C == *D);

    Interpolant* E = new Interpolant(max,
                                     xmin,
                                     xmax,
                                     max2,
                                     x2min,
                                     x2max,
                                     X_YY,
                                     romberg,
                                     rational,
                                     relative,
                                     isLog,
                                     romberg2,
                                     rational2,
                                     relative2,
                                     isLog2,
                                     rombergY,
                                     rationalY,
                                     relativeY,
                                     logSubst);
    Interpolant* F = new Interpolant(max,
                                     xmin,
                                     xmax,
                                     max2,
                                     x2min,
                                     x2max,
                                     X_YY,
                                     romberg,
                                     rational,
                                     relative,
                                     isLog,
                                     romberg2,
                                     rational2,
                                     relative2,
                                     isLog2,
                                     rombergY,
                                     rationalY,
                                     relativeY,
                                     logSubst);

    EXPECT_TRUE(*E == *F);

    SearchX        = 7;
    double SearchY = 11;

    PolValue = E->Interpolate(SearchX, SearchY);
    PolValue = F->Interpolate(SearchX, SearchY);

    EXPECT_TRUE(*E == *F);
}

TEST(Comparison, Comparison_not_equal)
{

    double SearchX = 5;
    double PolValue;
    Interpolant A;

    Interpolant* B = new Interpolant(
        max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);
    Interpolant* C = new Interpolant(
        max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);

    PolValue = C->Interpolate(SearchX);

    Interpolant* D = new Interpolant(max,
                                     xmin,
                                     xmax,
                                     max2,
                                     x2min,
                                     x2max,
                                     X_YY,
                                     romberg,
                                     rational,
                                     relative,
                                     isLog,
                                     romberg2,
                                     rational2,
                                     relative2,
                                     isLog2,
                                     rombergY,
                                     rationalY,
                                     relativeY,
                                     false);

    Interpolant* E = new Interpolant(max,
                                     xmin,
                                     xmax,
                                     max2,
                                     x2min,
                                     x2max,
                                     X_YY,
                                     romberg,
                                     rational,
                                     relative,
                                     isLog,
                                     romberg2,
                                     rational2,
                                     relative2,
                                     isLog2,
                                     rombergY,
                                     rationalY,
                                     relativeY,
                                     true);

    EXPECT_TRUE(A != *B);
    EXPECT_TRUE(*B != *C);
    EXPECT_TRUE(*D != *E);
}

TEST(Assignment, Copyconstructor)
{
    Interpolant A;
    Interpolant B = A;

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Copyconstructor2)
{
    Interpolant A(max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);
    Interpolant B(A);

    EXPECT_TRUE(A == B);

    Interpolant* C = new Interpolant(max,
                                     xmin,
                                     xmax,
                                     max2,
                                     x2min,
                                     x2max,
                                     X_YY,
                                     romberg,
                                     rational,
                                     relative,
                                     isLog,
                                     romberg2,
                                     rational2,
                                     relative2,
                                     isLog2,
                                     rombergY,
                                     rationalY,
                                     relativeY,
                                     true);
    Interpolant D(*C);

    EXPECT_TRUE(*C == D);
}

TEST(Assignment, Operator)
{
    Interpolant A(max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);
    Interpolant B;

    EXPECT_TRUE(A != B);

    B = A;

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Swap)
{
    Interpolant A(max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);
    Interpolant B;
    Interpolant* C = new Interpolant(max,
                                     xmin,
                                     xmax,
                                     max2,
                                     x2min,
                                     x2max,
                                     X_YY,
                                     romberg,
                                     rational,
                                     relative,
                                     isLog,
                                     romberg2,
                                     rational2,
                                     relative2,
                                     isLog2,
                                     rombergY,
                                     rationalY,
                                     relativeY,
                                     true);
    Interpolant D(max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);
    Interpolant E;
    Interpolant* F = new Interpolant(max,
                                     xmin,
                                     xmax,
                                     max2,
                                     x2min,
                                     x2max,
                                     X_YY,
                                     romberg,
                                     rational,
                                     relative,
                                     isLog,
                                     romberg2,
                                     rational2,
                                     relative2,
                                     isLog2,
                                     rombergY,
                                     rationalY,
                                     relativeY,
                                     true);
    EXPECT_TRUE(A == D);
    EXPECT_TRUE(B == E);
    EXPECT_TRUE(*C == *F);

    D.swap(E);

    EXPECT_TRUE(A == E);
    EXPECT_TRUE(B == D);

    D.swap(E);
    D.swap(*F);

    EXPECT_TRUE(*C == D);
}

TEST(_1D_Interpol, Simple_Test_of_XX_EXPX)
{
    Interpolant* Pol1 = new Interpolant(
        max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);
    double SearchX   = 5;
    double precision = 1E-5;

    double PolValue  = Pol1->Interpolate(SearchX);
    double RealValue = X2(SearchX);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    delete Pol1;
}

TEST(_1D_Interpol, Save_To_File)
{
    Interpolant* Pol1 = new Interpolant(
        max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);

    ASSERT_TRUE(Pol1->Save(File1DTest));

    delete Pol1;
}

TEST(_1D_Interpol, Load_From_File)
{
    Interpolant* Pol1 = new Interpolant();
    Pol1->Load(File1DTest);
    double SearchX   = 5;
    double precision = 1E-5;

    double PolValue  = Pol1->Interpolate(SearchX);
    double RealValue = X2(SearchX);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    delete Pol1;
}

TEST(_1D_Interpol, Rational_On)
{
    Interpolant* Pol1 =
        new Interpolant(max, xmin, xmax, X2, romberg, true, relative, isLog, rombergY, rationalY, relativeY, logSubst);
    double SearchX   = 5;
    double precision = 1E-5;

    double PolValue  = Pol1->Interpolate(SearchX);
    double RealValue = X2(SearchX);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    PolValue = Pol1->FindLimit(RealValue);
    ASSERT_NEAR(PolValue, SearchX, SearchX * precision);

    delete Pol1;
}

TEST(_1D_Interpol, isLog_On)
{
    Interpolant* Pol1 = new Interpolant(
        max, xmin, xmax, X2, romberg, rational, relative, true, rombergY, rationalY, relativeY, logSubst);
    double SearchX   = 5;
    double precision = 1E-5;

    double PolValue  = Pol1->Interpolate(SearchX);
    double RealValue = X2(SearchX);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    PolValue = Pol1->FindLimit(RealValue);
    ASSERT_NEAR(PolValue, SearchX, SearchX * precision);

    delete Pol1;
}

TEST(_1D_Interpol, RationalY_On)
{
    Interpolant* Pol1 = new Interpolant(
        max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, !rationalY, relativeY, logSubst);
    double SearchX   = 5;
    double precision = 1E-4;

    double PolValue  = Pol1->Interpolate(SearchX);
    double RealValue = X2(SearchX);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    PolValue = Pol1->FindLimit(RealValue);
    ASSERT_NEAR(PolValue, SearchX, SearchX * precision);

    delete Pol1;
}

TEST(_1D_Interpol, logSubst_On)
{
    Interpolant* Pol1 = new Interpolant(
        max, xmin, xmax, X2, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, !logSubst);
    double SearchX   = 5;
    double precision = 1E-5;

    double PolValue  = Pol1->Interpolate(SearchX);
    double RealValue = X2(SearchX);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    PolValue = Pol1->FindLimit(RealValue);
    ASSERT_NEAR(PolValue, SearchX, SearchX * precision);

    delete Pol1;
}

TEST(_1D_Interpol, All_On)
{
    Interpolant* Pol1 = new Interpolant(
        max, xmin, xmax, X2, romberg, !rational, relative, !isLog, rombergY, !rationalY, relativeY, !logSubst);
    double SearchX   = 5;
    double precision = 1E-5;

    double PolValue  = Pol1->Interpolate(SearchX);
    double RealValue = X2(SearchX);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    PolValue = Pol1->FindLimit(RealValue);
    ASSERT_NEAR(PolValue, SearchX, SearchX * precision);

    delete Pol1;
}

TEST(_2D_Interpol, Simple_Test_of_X_YY_EXPX)
{
    Interpolant* Pol2 = new Interpolant(max,
                                        xmin,
                                        xmax,
                                        max2,
                                        x2min,
                                        x2max,
                                        X_YY,
                                        romberg,
                                        rational,
                                        relative,
                                        isLog,
                                        romberg2,
                                        rational2,
                                        relative2,
                                        isLog2,
                                        rombergY,
                                        rationalY,
                                        relativeY,
                                        logSubst);
    double SearchX    = 7;
    double SearchY    = 11;
    double precision  = 1E-6;

    double PolValue  = Pol2->Interpolate(SearchX, SearchY);
    double RealValue = X_YY(SearchX, SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    delete Pol2;
}

TEST(_2D_Interpol, Save_To_File)
{
    Interpolant* Pol2 = new Interpolant(max,
                                        xmin,
                                        xmax,
                                        max2,
                                        x2min,
                                        x2max,
                                        X_YY,
                                        romberg,
                                        rational,
                                        relative,
                                        isLog,
                                        romberg2,
                                        rational2,
                                        relative2,
                                        isLog2,
                                        rombergY,
                                        rationalY,
                                        relativeY,
                                        logSubst);

    ASSERT_TRUE(Pol2->Save(File2DTest));

    delete Pol2;
}

TEST(_2D_Interpol, Load_From_File)
{
    Interpolant* Pol2 = new Interpolant();
    sleep(10);
    Pol2->Load(File2DTest);

    double SearchX   = 7;
    double SearchY   = 11;
    double precision = 1E-6;

    double PolValue  = Pol2->Interpolate(SearchX, SearchY);
    double RealValue = X_YY(SearchX, SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    delete Pol2;
}

TEST(_2D_Interpol, rational1_On)
{
    Interpolant* Pol2 = new Interpolant(max,
                                        xmin,
                                        xmax,
                                        max2,
                                        x2min,
                                        x2max,
                                        X_YY,
                                        romberg,
                                        !rational,
                                        relative,
                                        isLog,
                                        romberg2,
                                        rational2,
                                        relative2,
                                        isLog2,
                                        rombergY,
                                        rationalY,
                                        relativeY,
                                        logSubst);
    double SearchX    = 7;
    double SearchY    = 11;
    double precision  = 1E-6;

    double PolValue  = Pol2->Interpolate(SearchX, SearchY);
    double RealValue = X_YY(SearchX, SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    PolValue = Pol2->FindLimit(SearchX, RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY * precision);

    delete Pol2;
}

TEST(_2D_Interpol, rational2_On)
{
    Interpolant* Pol2 = new Interpolant(max,
                                        xmin,
                                        xmax,
                                        max2,
                                        x2min,
                                        x2max,
                                        X_YY,
                                        romberg,
                                        rational,
                                        relative,
                                        isLog,
                                        romberg2,
                                        !rational2,
                                        relative2,
                                        isLog2,
                                        rombergY,
                                        rationalY,
                                        relativeY,
                                        logSubst);
    double SearchX    = 7;
    double SearchY    = 11;
    double precision  = 1E-6;

    double PolValue  = Pol2->Interpolate(SearchX, SearchY);
    double RealValue = X_YY(SearchX, SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    PolValue = Pol2->FindLimit(SearchX, RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY * precision);

    delete Pol2;
}

TEST(_2D_Interpol, rational12_On)
{
    Interpolant* Pol2 = new Interpolant(max,
                                        xmin,
                                        xmax,
                                        max2,
                                        x2min,
                                        x2max,
                                        X_YY,
                                        romberg,
                                        !rational,
                                        relative,
                                        isLog,
                                        romberg2,
                                        !rational2,
                                        relative2,
                                        isLog2,
                                        rombergY,
                                        rationalY,
                                        relativeY,
                                        logSubst);
    double SearchX    = 7;
    double SearchY    = 11;
    double precision  = 1E-6;

    double PolValue  = Pol2->Interpolate(SearchX, SearchY);
    double RealValue = X_YY(SearchX, SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    PolValue = Pol2->FindLimit(SearchX, RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY * precision);

    delete Pol2;
}

TEST(_2D_Interpol, isLog1_On)
{
    Interpolant* Pol2 = new Interpolant(max,
                                        xmin,
                                        xmax,
                                        max2,
                                        x2min,
                                        x2max,
                                        X_YY,
                                        romberg,
                                        rational,
                                        relative,
                                        !isLog,
                                        romberg2,
                                        rational2,
                                        relative2,
                                        isLog2,
                                        rombergY,
                                        rationalY,
                                        relativeY,
                                        logSubst);
    double SearchX    = 7;
    double SearchY    = 11;
    double precision  = 1E-6;

    double PolValue  = Pol2->Interpolate(SearchX, SearchY);
    double RealValue = X_YY(SearchX, SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    PolValue = Pol2->FindLimit(SearchX, RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY * precision);

    delete Pol2;
}

TEST(_2D_Interpol, isLog2_On)
{
    Interpolant* Pol2 = new Interpolant(max,
                                        xmin,
                                        xmax,
                                        max2,
                                        x2min,
                                        x2max,
                                        X_YY,
                                        romberg,
                                        rational,
                                        relative,
                                        isLog,
                                        romberg2,
                                        rational2,
                                        relative2,
                                        !isLog2,
                                        rombergY,
                                        rationalY,
                                        relativeY,
                                        logSubst);
    double SearchX    = 7;
    double SearchY    = 11;
    double precision  = 1E-6;

    double PolValue  = Pol2->Interpolate(SearchX, SearchY);
    double RealValue = X_YY(SearchX, SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    PolValue = Pol2->FindLimit(SearchX, RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY * precision);

    delete Pol2;
}

TEST(_2D_Interpol, isLog12_On)
{
    Interpolant* Pol2 = new Interpolant(max,
                                        xmin,
                                        xmax,
                                        max2,
                                        x2min,
                                        x2max,
                                        X_YY,
                                        romberg,
                                        rational,
                                        relative,
                                        !isLog,
                                        romberg2,
                                        rational2,
                                        relative2,
                                        !isLog2,
                                        rombergY,
                                        rationalY,
                                        relativeY,
                                        logSubst);
    double SearchX    = 7;
    double SearchY    = 11;
    double precision  = 1E-6;

    double PolValue  = Pol2->Interpolate(SearchX, SearchY);
    double RealValue = X_YY(SearchX, SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    PolValue = Pol2->FindLimit(SearchX, RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY * precision);

    delete Pol2;
}

TEST(_2D_Interpol, rationalY_On)
{
    Interpolant* Pol2 = new Interpolant(max,
                                        xmin,
                                        xmax,
                                        max2,
                                        x2min,
                                        x2max,
                                        X_YY,
                                        romberg,
                                        rational,
                                        relative,
                                        isLog,
                                        romberg2,
                                        rational2,
                                        relative2,
                                        isLog2,
                                        rombergY,
                                        !rationalY,
                                        relativeY,
                                        logSubst);
    double SearchX    = 7;
    double SearchY    = 11;
    double precision  = 1E-6;

    double PolValue  = Pol2->Interpolate(SearchX, SearchY);
    double RealValue = X_YY(SearchX, SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    PolValue = Pol2->FindLimit(SearchX, RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY * precision);

    delete Pol2;
}

TEST(_2D_Interpol, logSubst_On)
{
    Interpolant* Pol2 = new Interpolant(max,
                                        xmin,
                                        xmax,
                                        max2,
                                        x2min,
                                        x2max,
                                        X_YY,
                                        romberg,
                                        rational,
                                        relative,
                                        isLog,
                                        romberg2,
                                        rational2,
                                        relative2,
                                        isLog2,
                                        rombergY,
                                        rationalY,
                                        relativeY,
                                        !logSubst);
    double SearchX    = 7;
    double SearchY    = 11;
    double precision  = 1E-6;

    double PolValue  = Pol2->Interpolate(SearchX, SearchY);
    double RealValue = X_YY(SearchX, SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    PolValue = Pol2->FindLimit(SearchX, RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY * precision);

    delete Pol2;
}

TEST(_2D_Interpol, All_On)
{
    Interpolant* Pol2 = new Interpolant(max,
                                        xmin,
                                        xmax,
                                        max2,
                                        x2min,
                                        x2max,
                                        X_YY,
                                        romberg,
                                        !rational,
                                        relative,
                                        !isLog,
                                        romberg2,
                                        !rational2,
                                        relative2,
                                        !isLog2,
                                        rombergY,
                                        !rationalY,
                                        relativeY,
                                        !logSubst);
    double SearchX    = 7;
    double SearchY    = 11;
    double precision  = 1E-6;

    double PolValue  = Pol2->Interpolate(SearchX, SearchY);
    double RealValue = X_YY(SearchX, SearchY);

    ASSERT_NEAR(PolValue, RealValue, RealValue * precision);

    PolValue = Pol2->FindLimit(SearchX, RealValue);
    ASSERT_NEAR(PolValue, SearchY, SearchY * precision);

    delete Pol2;
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();

    // Clean the mess
    // std::remove (File1DTest.c_str());
    // std::remove (File2DTest.c_str());
    return ret;
}
