#include "gtest/gtest.h"

#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/Time.h"
#include "PROPOSAL/crossection/CrossSectionBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"

using namespace PROPOSAL;

TEST(ApproximateTimeBuilder, Constructor) {
    Time* time_approx1 = new ApproximateTimeBuilder();
    ApproximateTimeBuilder time_approx2();
}

TEST(ExactTimeBuilder, Constructor) {
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy", MuMinusDef());
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});

    Time* time_exact1 = new ExactTimeBuilder<UtilityIntegral>(CrossSectionList{cross_dummy}, MuMinusDef());
    Time* time_exact2 = new ExactTimeBuilder<UtilityInterpolant>(CrossSectionList{cross_dummy}, MuMinusDef());

    ExactTimeBuilder<UtilityIntegral> time_exact3(CrossSectionList{cross_dummy}, MuMinusDef());
    ExactTimeBuilder<UtilityInterpolant> time_exact4(CrossSectionList{cross_dummy}, MuMinusDef());
}

class ExactTimeBuilderDummy : public ExactTimeBuilder<UtilityIntegral>{
public:
    ExactTimeBuilderDummy(CrossSectionList cross, const ParticleDef& p_def) : ExactTimeBuilder<UtilityIntegral>(cross, p_def){};
    double GetMass(){return mass;}
};

TEST(ExactTimeBuilder, ConstructorMass){
    auto cross_dummy_MMU = std::make_shared<CrossSectionBuilder>("contrand_dummy", MuMinusDef());
    cross_dummy_MMU->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});
    auto time_exact_MMU = ExactTimeBuilderDummy(CrossSectionList{cross_dummy_MMU}, MuMinusDef());
    EXPECT_DOUBLE_EQ(time_exact_MMU.GetMass(), MMU);

    auto cross_dummy_ME = std::make_shared<CrossSectionBuilder>("contrand_dummy", EMinusDef());
    cross_dummy_ME->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});
    auto time_exact_ME = ExactTimeBuilderDummy(CrossSectionList{cross_dummy_ME}, EMinusDef());
    EXPECT_DOUBLE_EQ(time_exact_ME.GetMass(), ME);
}

TEST(ApproximateTimeBuilder, ZeroDistance){
    Time* time_approx1 = new ApproximateTimeBuilder();
    EXPECT_DOUBLE_EQ(time_approx1->TimeElapsed(0), 0.);
}

TEST(ExactTimeBuilder, EqualEnergies){
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy", MuMinusDef());
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});
    auto time_exact = ExactTimeBuilder<UtilityIntegral>(CrossSectionList{cross_dummy}, MuMinusDef());

    EXPECT_DOUBLE_EQ(time_exact.TimeElapsed(1e5, 1e5, 0.), 0.);
}

TEST(ApproximateTimeBuilder, PhysicalBehaviour){
    auto time_approx = ApproximateTimeBuilder();
    auto distances = std::array<double, 4>{1e2, 1e4, 1e6, 1e8};

    double elapsed_time_old = 0;
    double aux;
    for(auto d : distances){
        aux = time_approx.TimeElapsed(d);
        EXPECT_TRUE(elapsed_time_old < aux);
        elapsed_time_old = aux;
    }
}

TEST(ExactTimeBuilder, PhysicalBehaviour){
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy", MuMinusDef());
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});
    auto time_exact = ExactTimeBuilder<UtilityIntegral>(CrossSectionList{cross_dummy}, MuMinusDef());

    double E_f = 1e3;
    auto energies = std::array<double, 9>{1e4, 1e5, 1e6, 1e7,
                                          1e8, 1e9, 1e10, 1e11, 1e12};

    double elapsed_time_old = 0;
    double aux;
    for(auto E_i : energies){
        aux = time_exact.TimeElapsed(E_i, E_f, 0.);
        EXPECT_TRUE(elapsed_time_old < aux);
        elapsed_time_old = aux;
    }
}

TEST(TimeBuilder, HighEnergies){
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy", MuMinusDef());
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});

    auto time_exact = ExactTimeBuilder<UtilityIntegral>(CrossSectionList{cross_dummy}, MuMinusDef());
    auto time_approx = ApproximateTimeBuilder();

    auto displacement = DisplacementBuilder<UtilityIntegral>(CrossSectionList{cross_dummy});

    double elapsed_time_exact = time_exact.TimeElapsed(1e13, 1e12, 0);
    double prop_distance = displacement.SolveTrackIntegral(1e13, 1e12, 0);

    double elapsed_time_approx = time_approx.TimeElapsed(prop_distance);

    EXPECT_FALSE(elapsed_time_approx == elapsed_time_exact);
    EXPECT_TRUE(elapsed_time_approx < elapsed_time_exact); //exact velocity must be smaller than c
    EXPECT_NEAR(elapsed_time_approx, elapsed_time_exact, elapsed_time_exact*1e-5);
}

TEST(TimeBuilder, MasslessParticles){
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy", GammaDef());
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});

    auto time_exact = ExactTimeBuilder<UtilityIntegral>(CrossSectionList{cross_dummy}, GammaDef());
    auto time_approx = ApproximateTimeBuilder();

    auto displacement = DisplacementBuilder<UtilityIntegral>(CrossSectionList{cross_dummy});

    double elapsed_time_exact = time_exact.TimeElapsed(1e6, 1e4, 0);
    double prop_distance = displacement.SolveTrackIntegral(1e6, 1e4, 0);

    double elapsed_time_approx = time_approx.TimeElapsed(prop_distance);

    EXPECT_NEAR(elapsed_time_approx,  elapsed_time_exact, elapsed_time_exact*1e-10); //massless particles should have a velocity of c
}

TEST(ExactTimeBuilder, compare_integral_interpolant){
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy", MuMinusDef());
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});

    auto time_integral = ExactTimeBuilder<UtilityIntegral>(CrossSectionList{cross_dummy}, MuMinusDef());
    auto time_interpolant = ExactTimeBuilder<UtilityInterpolant>(CrossSectionList{cross_dummy}, MuMinusDef());

    double E_f = 1e5;
    auto energies = std::array<double, 4>{1e7, 1e9, 1e11, 1e13};
    double integrated_time, interpolated_time;
    for(auto E_i : energies){
        integrated_time = time_integral.TimeElapsed(E_i, E_f, 0.);
        interpolated_time = time_interpolant.TimeElapsed(E_i, E_f, 0.);
        EXPECT_NEAR(integrated_time, interpolated_time, integrated_time * 1e-5);
        EXPECT_FALSE(integrated_time == interpolated_time);
    }
}

TEST(ExactTimeBuilder, TimeElapsedUpperLimit){
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy", MuMinusDef());
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});
    auto exact_time = ExactTimeBuilder<UtilityIntegral>(CrossSectionList{cross_dummy}, MuMinusDef());

    double E_f = 1e4;
    auto energies = std::array<double, 4>{1e6, 1e8, 1e10, 1e12};

    for(auto E_i : energies){
        double elapsed_time = exact_time.TimeElapsed(E_i, E_f, 0);
        double upper_limit = exact_time.TimeElapsedUpperLimit(E_i, elapsed_time);
        EXPECT_NEAR(upper_limit, E_f, E_f * 1e-4);
    }
}

TEST(ExactTimeBuilder, TimeElapsedUpperLimit_interpolant){
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy", MuMinusDef());
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});

    auto time_integral = ExactTimeBuilder<UtilityIntegral>(CrossSectionList{cross_dummy}, MuMinusDef());
    auto time_interpol = ExactTimeBuilder<UtilityInterpolant>(CrossSectionList{cross_dummy}, MuMinusDef());

    auto energies = std::array<double, 4>{1e6, 1e8, 1e10, 1e12};

    for(auto E_i : energies){
        double upper_limit_integral = time_integral.TimeElapsedUpperLimit(E_i, 1e-7);
        double upper_limit_interpol = time_interpol.TimeElapsedUpperLimit(E_i, 1e-7);
        EXPECT_NEAR(upper_limit_integral, upper_limit_interpol, upper_limit_integral * 1e-5);
    }
}

TEST(ExactTimeBuilder, TimeElapsedUpperLimit_zero){
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy", MuMinusDef());
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});
    auto exact_time = ExactTimeBuilder<UtilityIntegral>(CrossSectionList{cross_dummy}, MuMinusDef());

    EXPECT_DOUBLE_EQ(exact_time.TimeElapsedUpperLimit(1e6, 0.), 1e6);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
