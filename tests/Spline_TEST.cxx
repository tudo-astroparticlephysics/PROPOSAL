#include <cmath>
#include <iostream>
#include <numeric>

#include "gtest/gtest.h"

#include "PROPOSAL/PROPOSAL.h"

using namespace PROPOSAL;

double cubic(double x) {
    double result = 0.0;
    for (int i = 0; i < 3; ++i) {
        result += std::pow(x, i);
    }
    return result;
}

double linear(double x) {
    return x + 1;
}


TEST(Linear_Spline, Comparison_equal) {
    std::vector<double> x(10);
    std::vector<double> f_x(10);
    std::iota(std::begin(x), std::end(x), 0);

    for (int i = 0; i < x.size(); ++i) {
        f_x[i] = linear(x[i]);
    }
    Linear_Spline sp(x, f_x);
    for (int i = 0; i < x.size(); ++i) {
        ASSERT_NEAR(sp.evaluate(x[i]), f_x[i], COMPUTER_PRECISION * f_x[i]);
    }
}

TEST(Cubic_Spline, Comparison_equal) {
    std::vector<double> x(10);
    std::vector<double> f_x(10);
    std::iota(std::begin(x), std::end(x), 0);

    for (int i = 0; i < x.size(); ++i) {
        f_x[i] = cubic(x[i]);
    }
    Cubic_Spline sp(x, f_x);
    for (int i = 0; i < x.size(); ++i) {
        ASSERT_NEAR(sp.evaluate(x[i]), f_x[i], COMPUTER_PRECISION * f_x[i]);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
