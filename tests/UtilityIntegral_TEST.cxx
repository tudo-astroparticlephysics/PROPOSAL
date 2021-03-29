#include "gtest/gtest.h"
#include <cmath>
#include <array>

#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"

using namespace PROPOSAL;

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

constexpr static size_t hash = 42;
TEST(Constructor, Constructor){
    auto integrand = [](double x)->double {return 1/x;};
    auto integral = UtilityIntegral(integrand, 0, hash);
}

TEST(Calculate, Integration){
    auto integrand = [](double x)->double {return 1/x;};
    auto integral = UtilityIntegral(integrand, 0, hash);

    for(double logmin = -7; logmin <= 7; logmin+=1e-1){
        for(double logmax = -7; logmax <= 7; logmax+=1e-1){
            double min = std::pow(10., logmin);
            double max = std::pow(10., logmax);
            auto analytical = std::log(max) - std::log(min);
            EXPECT_NEAR(integral.Calculate(min, max), analytical, std::abs(analytical * 1e-5));
        }
    }
}

TEST(Calculate, BackwardIntegration){
    auto integrand = [](double x)->double {return -1/x;};
    auto integral = UtilityIntegral(integrand, 0, hash);

    for(double logmin = -7; logmin <= 7; logmin+=1e-1){
        for(double logmax = -7; logmax <= 7; logmax+=1e-1){
            double min = std::pow(10., logmin);
            double max = std::pow(10., logmax);
            auto analytical = std::log(min) - std::log(max);
            EXPECT_NEAR(integral.Calculate(min, max), analytical, std::abs(analytical * 1e-5));
        }
    }
}

TEST(GetUpperLimit, ForwardIntegration){
    //Test forward integration with positive integrand
    auto integrand1 = [](double x)->double { return 1/x;};
    double HIGH = 1.e14;
    auto integral1 = UtilityIntegral(integrand1, HIGH, hash);

    for (double logmin = -7; logmin < 7; logmin+=1e-1) {
        for (double logxi = -3; logxi < 3; logxi+=1e-1) {
            double min = std::pow(10., logmin);
            double xi = std::pow(10., logxi);
            double analytical_max = min * std::exp(xi);
            if (analytical_max > HIGH) {
                #ifndef NDEBUG
                EXPECT_DEATH(integral1.GetUpperLimit(min, xi), "");
                #endif
            } else {
                EXPECT_NEAR(integral1.GetUpperLimit(min, xi), analytical_max, analytical_max*1e-5);
            }
        }
    }
}

TEST(GetUpperLimit, ForwardIntegration_Outside_Range){
    //Test GetUpperLimit when expected upper_limit is bigger than high parameter
    auto integrand1 = [](double x)->double { return 1/x;};
    double high = 1e14;
    auto integral1 = UtilityIntegral(integrand1, high, hash);

    double lower_lim = 1e2;
    double xi = 40;

    EXPECT_TRUE(lower_lim * std::exp(xi) > high);
    #ifndef NDEBUG
    EXPECT_DEATH(integral1.GetUpperLimit(lower_lim, xi), "");
    #endif
}

TEST(GetUpperLimit, BackwardIntegration){
    //Test backward integration with negative integrand
    auto integrand1 = [](double x)->double { return -1/x;};
    double LOW = 1;
    auto integral1 = UtilityIntegral(integrand1, LOW, hash);

    for (double loglower = 1; loglower < 7; loglower+=1e-1) {
        for (double logxi = -3; logxi < 3; logxi+=1e-1) {
            double lower = std::pow(10, loglower);
            double xi = std::pow(10, logxi);
            double analytical_upper = lower * std::exp(-xi);
            if (analytical_upper < LOW) {
                #ifndef NDEBUG
                EXPECT_DEATH(integral1.GetUpperLimit(lower, xi), "");
                #endif
            } else {
                EXPECT_NEAR(integral1.GetUpperLimit(lower, xi), analytical_upper, analytical_upper*1e-5);
            }
        }
    }
}

TEST(GetUpperLimit, BackwardIntegration_Outside_Range){
    //Test GetUpperLimit when expected lower_lim is smaller than low parameter
    auto integrand1 = [](double x)->double { return -1/x;};
    double LOW = 10;
    auto integral1 = UtilityIntegral(integrand1, LOW, hash);
    double xi = 10;

    EXPECT_TRUE(LOW * std::exp(-xi) < LOW);
    #ifndef NDEBUG
    double lower = 1e5;
    EXPECT_DEATH(integral1.GetUpperLimit(lower, xi), "");
    #endif
}

TEST(GetUpperLimit, ListForwardIntegration){
    auto integrand1 = [](double x)->double { return 1/x;};
    auto integral = UtilityIntegral(integrand1, 1e14, hash);

    for(double loglower = 1.; loglower < 13.; loglower+=1e-1){
        for(auto logupper = 1.; logupper < 13.; logupper+=1e-1){
            double lower = std::pow(10, loglower);
            double upper = std::pow(10, logupper);
            if(lower > upper){
                continue;
            }
            double xi = integral.Calculate(lower, upper); //Calculate result
            EXPECT_NEAR(integral.GetUpperLimit(lower, xi), upper, upper * 1e-5);
        }
    }
}

TEST(GetUpperLimit, ListBackwardIntegration){
    auto integrand1 = [](double x)->double { return -1/x;};
    auto integral = UtilityIntegral(integrand1, 1, hash);

    for(double loglower = 1.; loglower < 13; loglower+=1e-1){
        for(double logupper = 1; logupper < 13; logupper+=1e-1){
            double lower = std::pow(10, loglower);
            double upper = std::pow(10, logupper);
                if(lower < upper){
                    continue;
                }
                double xi = integral.Calculate(lower, upper); //Calculate result
                EXPECT_NEAR(integral.GetUpperLimit(lower, xi), upper, upper * 1e-5);
       }
    }
}
