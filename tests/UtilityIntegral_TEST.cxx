#include "gtest/gtest.h"

#include "PROPOSAL/propagation_utility/PropagationUtility.h"

using namespace PROPOSAL;

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(Constructor, Constructor){
    auto integrand = [](double x)->double {return 1/x;};
    auto integral = UtilityIntegral(integrand, 0);
}

TEST(Calculate, Integration){
    auto integrand = [](double x)->double {return 1/x;};
    auto integral = UtilityIntegral(integrand, 0);
    auto limits = std::array<double, 7>{1e-7, 1e-3, 1, 10, 1e3, 1e5, 1e10};

    for(auto min : limits){
        for(auto max : limits){
            auto analytical = std::log(max) - std::log(min);
            EXPECT_NEAR(integral.Calculate(min, max, 0.), analytical, std::abs(analytical * 1e-5));
        }
    }
}

TEST(GetUpperLimit, ForwardIntegration){
    //Test forward integration with positive integrand
    auto integrand1 = [](double x)->double { return 1/x;};
    auto integral1 = UtilityIntegral(integrand1, 1e14);

    double lower_lim = 1;
    double xi = 10;
    double analytical_upper = lower_lim * std::exp(xi);

    EXPECT_NEAR(integral1.GetUpperLimit(lower_lim, xi), analytical_upper, analytical_upper*1e-5);
}

TEST(GetUpperLimit, ForwardIntegration_Outside_Range){
    //Test GetUpperLimit when expected upper_limit is bigger than high parameter (expect high to be returned)
    auto integrand1 = [](double x)->double { return 1/x;};
    double high = 1e14;
    auto integral1 = UtilityIntegral(integrand1, high);

    double lower_lim = 1e2;
    double xi = 40;

    EXPECT_NEAR(integral1.GetUpperLimit(lower_lim, xi), high, high*1e-5);
}

TEST(GetUpperLimit, BackwardIntegration){
    //Test backward integration with negative integrand
    auto integrand1 = [](double x)->double { return -1/x;};
    auto integral1 = UtilityIntegral(integrand1, 1);

    double lower_lim = 1e6;
    double xi = 5;
    double analytical_upper = lower_lim * std::exp(-xi);

    EXPECT_NEAR(integral1.GetUpperLimit(lower_lim, xi), analytical_upper, analytical_upper*1e-5);
}

TEST(GetUpperLimit, BackwardIntegration_Outside_Range){
    //Test GetUpperLimit when expected lower_lim is smaller than low parameter (expect low to be returned)
    auto integrand1 = [](double x)->double { return -1/x;};
    double low = 10;
    auto integral1 = UtilityIntegral(integrand1, low);

    double upper_lim = 1e5;
    double xi = 10;

    EXPECT_NEAR(integral1.GetUpperLimit(upper_lim, xi), low, low*1e-5);
}

TEST(GetUpperLimit, ListForwardIntegration){
    auto integrand1 = [](double x)->double { return 1/x;};
    auto integral = UtilityIntegral(integrand1, 1e14);
    auto limits = std::array<double, 6>{1e3, 1e5, 1e7, 1e9, 1e11, 1e13};

    for(auto lower_lim : limits){
        for(auto upper_lim : limits){
            if(lower_lim >= upper_lim){
                continue;
            }
            double xi = integral.Calculate(lower_lim,  upper_lim, 0.); //Calculate result
            EXPECT_NEAR(integral.GetUpperLimit(lower_lim, xi), upper_lim, upper_lim * 1e-5);
        }
    }
}

TEST(GetUpperLimit, ListBackwardIntegration){
    auto integrand1 = [](double x)->double { return -1/x;};
    auto integral = UtilityIntegral(integrand1, 1);
    auto limits = std::array<double, 4>{1e3, 1e4, 1e5, 1e6};

    for(auto E_i : limits){
        for(auto E_f : limits){
                if(E_i <= E_f){
                    continue;
                }
                double xi = integral.Calculate(E_i,  E_f, 0.); //Calculate result
                EXPECT_NEAR(integral.GetUpperLimit(E_i, xi), E_f, E_f * 1e-5);
       }
    }
}