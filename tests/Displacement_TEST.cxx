#include "gtest/gtest.h"

#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include "PROPOSAL/crosssection/ParticleDefaultCrossSectionList.h"
#include <vector>

using std::vector;
using namespace PROPOSAL;

auto GetCrossSections() {
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 1, false);
    auto cross = GetStdCrossSections(MuMinusDef(), Ice(), cuts, true);
    return cross;
}

TEST(constructor, integral)
{
    DisplacementBuilder disp(GetCrossSections(), std::false_type());
}

TEST(constructor, interpolant)
{
    DisplacementBuilder disp(GetCrossSections(), std::true_type());
}

TEST(constructor, error)
{
    auto cross = crosssection_list_t();
    ASSERT_THROW(DisplacementBuilder disp(cross, std::false_type()),
        std::invalid_argument);
}

TEST(SolveTrackIntegral, ConsistencyCheck)
{
    // If the energy difference increases, the result of SolveTrack integral
    // (e.g. the grammage) should increase as well
    DisplacementBuilder disp(GetCrossSections(), std::false_type());
    double E_f = 1e3;
    double displacement;
    double displacement_previous = -1.;
    for (double logE_i = 3.; logE_i < 8; logE_i+=1.e-3) {
        displacement = disp.SolveTrackIntegral(std::pow(logE_i, 10.), E_f);
        EXPECT_GT(displacement, displacement_previous);
        displacement_previous = displacement;
    }
}

TEST(SolveTrackIntegral, ZeroDisplacement)
{
    // no energy difference should equal to no displacement
    DisplacementBuilder disp(GetCrossSections(), std::false_type());
    EXPECT_TRUE(disp.SolveTrackIntegral(1e6, 1e6) == 0.);
}

TEST(SolveTrackIntegral, CompareIntegralInterpolant)
{
    // comparing intergral and interpolant values
    DisplacementBuilder disp_calc_integral(GetCrossSections(), std::false_type());
    DisplacementBuilder disp_calc_interpol(GetCrossSections(), std::true_type());
    double E_f = 1e3;
    for (double logE_i = 3.; logE_i < 8; logE_i+=1.e-3) {
        double E_i = std::pow(logE_i, 10.);
        auto disp_integral = disp_calc_integral.SolveTrackIntegral(E_i, E_f);
        auto disp_interpol = disp_calc_interpol.SolveTrackIntegral(E_i, E_f);
        EXPECT_NEAR(disp_integral, disp_interpol, disp_integral*1e-3);
    }
}

TEST(UpperLimitTrackIntegral, ConsistencyCheck)
{
    // The final energy should become smaller for increasing displacements
    DisplacementBuilder disp_calc(GetCrossSections(), std::false_type());
    double E_i = 1.e8;
    double E_f_old = E_i;
    for (double logDisp = 1; logDisp < 5; logDisp+=1.e-2) {
        double disp = std::pow(logDisp, 10.);
        double E_f = disp_calc.UpperLimitTrackIntegral(E_i, disp);
        EXPECT_LE(E_f, E_f_old);
        E_f_old = E_f;
    }
}

TEST(UpperLimitTrackIntegral, ZeroDisplacement)
{
    DisplacementBuilder disp_calc(GetCrossSections(), std::false_type());
    EXPECT_NEAR( disp_calc.UpperLimitTrackIntegral(1e6, 0.), 1e6, 1e-6);
}

TEST(UpperLimitTrackIntegral, CompareIntegralInterpolant)
{
    // comparing intergral and interpolant values
    DisplacementBuilder disp_calc_integral(GetCrossSections(), std::false_type());
    DisplacementBuilder disp_calc_interpol(GetCrossSections(), std::true_type());
    double E_i = 1e8;
    for (double logDisp = 1; logDisp < 5; logDisp+=1.e-2) {
        double disp = std::pow(logDisp, 10.);
        double E_f_integral = disp_calc_integral.UpperLimitTrackIntegral(E_i, disp);
        double E_f_interpol;
        try {
            E_f_interpol = disp_calc_interpol.UpperLimitTrackIntegral(E_i, disp);
        } catch (std::logic_error& e) {
            E_f_interpol = disp_calc_interpol.GetLowerLim();
        }
        EXPECT_NEAR(E_f_integral, E_f_interpol, E_f_integral*1e-3);
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
