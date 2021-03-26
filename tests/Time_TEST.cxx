#include "gtest/gtest.h"

#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/crosssection/ParticleDefaultCrossSectionList.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include "PROPOSAL/propagation_utility/TimeBuilder.h"

using namespace PROPOSAL;

auto GetCrossSections(double interpol)
{
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 1, false);
    auto cross = GetStdCrossSections(MuMinusDef(), Ice(), cuts, interpol);
    return cross;
}

auto GetCrossSectionsGamma(double interpol)
{
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 1, false);
    auto cross = GetStdCrossSections(GammaDef(), Ice(), cuts, interpol);
    return cross;
}

TEST(ApproximateTimeBuilder, Constructor)
{
    Time* time_approx1 = new ApproximateTimeBuilder();
    ApproximateTimeBuilder time_approx2;

    delete time_approx1;
}

TEST(ExactTimeBuilder, Constructor)
{
    // constructor muon
    auto cross_mu_integral = GetCrossSections(false);
    auto cross_mu_interpol = GetCrossSections(true);

    auto time_exact1 = make_time(cross_mu_integral, MuMinusDef(), false);
    auto time_exact2 = make_time(cross_mu_interpol, MuMinusDef(), true);

    // constructor gamma
    auto cross_gamma_integral = GetCrossSectionsGamma(false);
    auto cross_gamma_interpol = GetCrossSectionsGamma(true);

    auto time_exact3 = make_time(cross_gamma_integral, GammaDef(), false);
    auto time_exact4 = make_time(cross_gamma_interpol, GammaDef(), true);
}

TEST(ApproximateTimeBuilder, ZeroDistance)
{
    // No time should elapse if no distance is propagated
    ApproximateTimeBuilder time_approx;
    EXPECT_DOUBLE_EQ(time_approx.TimeElapsed(1e6, 1e6, 0., 1.), 0.);
}

TEST(ExactTimeBuilder, ZeroDistance)
{
    // No time should elapse if no distance is propagated
    auto cross_integral = GetCrossSections(false);
    auto cross_interpol = GetCrossSections(true);

    auto exact_time1 = make_time(cross_integral, MuMinusDef(), false);
    auto exact_time2 = make_time(cross_interpol, MuMinusDef(), true);
    EXPECT_DOUBLE_EQ(exact_time1->TimeElapsed(1e6, 1e6, 0., 1.), 0.);
    EXPECT_DOUBLE_EQ(exact_time2->TimeElapsed(1e6, 1e6, 0., 1.), 0.);
}

TEST(ApproximateTimeBuilder, PhysicalBehaviour)
{
    // with increasing distance, the elapsed time should increase as well

    auto time_approx = ApproximateTimeBuilder();
    double elapsed_time_old = 0;
    for (double log_dist = 1; log_dist < 8; log_dist += 0.1) {
        double dist = std::pow(log_dist, 10.);
        double elapsed_time = time_approx.TimeElapsed(0., 0., dist, 1.);
        EXPECT_GT(elapsed_time, elapsed_time_old);
        elapsed_time_old = elapsed_time;
    }
}

TEST(ExactTimeBuilder, PhysicalBehaviour)
{
    // with increasing energy difference, the elapsed time should increase as
    // well

    auto cross = GetCrossSections(true);
    auto medium = Ice();
    auto time_exact = make_time(cross, MuMinusDef(), true);
    double E_f = 1e3;

    double elapsed_time_old = -1;
    for (double Elog_i = 3.; Elog_i < 12.; Elog_i += 0.1) {
        double E_i = std::pow(Elog_i, 10.);
        double elapsed_time
            = time_exact->TimeElapsed(E_i, E_f, 0., medium.GetMassDensity());
        EXPECT_GT(elapsed_time, elapsed_time_old);
        elapsed_time_old = elapsed_time;
    }
}

TEST(TimeBuilder, HighEnergies)
{
    auto cross = GetCrossSections(false);
    auto medium = Ice();

    auto time_exact = make_time(cross, MuMinusDef(), false);
    auto time_approx = ApproximateTimeBuilder();

    auto displacement = make_displacement(cross, false);

    double elapsed_time_exact
        = time_exact->TimeElapsed(1e14, 1e12, 0, medium.GetMassDensity());
    double prop_grammage = displacement->SolveTrackIntegral(1e14, 1e12);

    double elapsed_time_approx = time_approx.TimeElapsed(
        1e14, 1e12, prop_grammage, medium.GetMassDensity());

    EXPECT_LT(elapsed_time_approx, elapsed_time_exact); // exact velocity must be smaller than c
    EXPECT_NE(elapsed_time_approx, elapsed_time_exact); // they should not be equal
    EXPECT_NEAR(elapsed_time_approx, elapsed_time_exact,
        elapsed_time_exact
            * 1e-5); //... but they should be pretty close to equal
}

TEST(TimeBuilder, MasslessParticles)
{
    auto cross = GetCrossSectionsGamma(false);
    auto medium = Ice();

    auto time_exact = make_time(cross, GammaDef(), false);
    auto time_approx = ApproximateTimeBuilder();

    auto displacement = DisplacementBuilder(cross, std::false_type());
    double prop_grammage = displacement.SolveTrackIntegral(1e6, 1e4);

    double elapsed_time_exact = time_exact->TimeElapsed(
        1e6, 1e4, prop_grammage, medium.GetMassDensity());

    double elapsed_time_approx = time_approx.TimeElapsed(
        1e6, 1e4, prop_grammage, medium.GetMassDensity());

    EXPECT_NEAR(elapsed_time_approx, elapsed_time_exact,
        elapsed_time_exact
            * 1e-10); // massless particles should have a velocity of c
}

TEST(ExactTimeBuilder, CompareIntegralInterpolant)
{
    auto cross = GetCrossSections(true);
    auto medium = Ice();

    auto time_integral = make_time(cross, MuMinusDef(), false);
    auto time_interpolant = make_time(cross, MuMinusDef(), true);

    double E_f = 1e5;
    double integrated_time, interpolated_time;
    for (double Elog_i = 5.; Elog_i < 10; Elog_i += 1.e-2) {
        double E_i = std::pow(Elog_i, 10.);
        integrated_time
            = time_integral->TimeElapsed(E_i, E_f, 0., medium.GetMassDensity());
        interpolated_time = time_interpolant->TimeElapsed(
            E_i, E_f, 0., medium.GetMassDensity());
        EXPECT_NEAR(integrated_time, interpolated_time, integrated_time * 1e-5);
        EXPECT_FALSE(integrated_time == interpolated_time);
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
