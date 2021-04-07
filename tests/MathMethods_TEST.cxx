
#include <cmath>
#include <iostream>

#include "gtest/gtest.h"

#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/RandomGenerator.h"

using namespace PROPOSAL;

double f1(double x){
    return pow(x, 3) - 4;
}

double df1(double x){
    return 3. * pow(x, 2);
}

double f2(double x){
    return std::sin(x);
}

double df2(double x){
    return std::cos(x);
}


TEST(NewtonRaphson, Comparison_equal)
{

double real1 = pow(2., 2./3);
double real2 = 0;

double aux1 = NewtonRaphson(f1, df1, -5, 15, 5);
double aux2 = NewtonRaphson(f2, df2, -3./4 * PI, 1./2 * PI, -1./8 * PI);

ASSERT_NEAR(aux1, real1, 1.e-6 * real1);
ASSERT_NEAR(aux2, 0, 1.e-6 * aux2);

aux1 = NewtonRaphson(f1, df1, -5, 15, 15);
aux2 = NewtonRaphson(f2, df2, -3./4 * PI, 1./2 * PI, 1./2 * PI);

ASSERT_NEAR(aux1, real1, 1.e-6 * real1);
ASSERT_NEAR(aux2, real2, 1.e-6 * aux2);

}

TEST(Bisection, Comparison_equal)
{

    double real1 = pow(2., 2./3);
    double real2 = 0;

    double aux1 = Bisection(f1, -5, 15, 1e-6, 100);
    double aux2 = Bisection(f2, -3./4 * PI, 1./2 * PI, 1e-6, 100);

    ASSERT_NEAR(aux1, real1, 1.e-6 );
    ASSERT_NEAR(aux2, 0, 1.e-6 );

    aux1 = Bisection(f1, 15, -5, 1e-6, 100);
    aux2 = Bisection(f2, 1./2 * PI, -3./4 * PI, 1e-6, 100);

    ASSERT_NEAR(aux1, real1, 1.e-6 );
    ASSERT_NEAR(aux2, real2, 1.e-6 );

}


TEST(SampleFromGaussian, Momenta){
    RandomGenerator::Get().SetSeed(24601);
    double sigma = 2;
    double mean = 5;
    unsigned int statistics = 1e6;
    auto average = std::pair<double, double>{0., 0.};

    for(unsigned int n=1; n<=statistics; n++){
        double sampled = SampleFromGaussian(mean, sigma, RandomGenerator::Get().RandomDouble());
        average = welfords_online_algorithm(sampled, n, average.first, average.second);
    }
    EXPECT_NEAR(average.first, mean, 1e-3 * mean);
    EXPECT_NEAR(std::sqrt(average.second), sigma, 1e-3 * sigma);

}

TEST(SampleFromGaussian, Inf){
    RandomGenerator::Get().SetSeed(24601);
    double sigma = 0;
    double mean = 5;
    double sampled = SampleFromGaussian(mean, sigma, RandomGenerator::Get().RandomDouble());
    EXPECT_DOUBLE_EQ(mean, sampled);
}

TEST(SampleFromExponential, Momenta){
    RandomGenerator::Get().SetSeed(24601);
    double lambda = 10;
    double mean = 1./lambda;
    double sigma = 1./lambda;
    int statistics = 1e7;
    auto average = std::pair<double, double>{0., 0.};

    for(unsigned int n=1; n<=statistics; n++){
        double sampled = SampleFromExponential(RandomGenerator::Get().RandomDouble(), lambda);
        average = welfords_online_algorithm(sampled, n, average.first, average.second);
    }
    EXPECT_NEAR(average.first, mean, 1e-3 * mean);
    EXPECT_NEAR(std::sqrt(average.second), sigma, 1e-3 * sigma);

}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
