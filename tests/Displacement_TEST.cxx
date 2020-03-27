#include "gtest/gtest.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/propagation_utility/Displacement.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include <vector>

using std::vector;
using namespace PROPOSAL;

// dEdx = a + b * E
constexpr auto a = 1.e-1;
constexpr auto b = 1.e-3;

TEST(constructor, integral)
{
    CrossSectionList crosss{ std::make_shared<CrossSectionBuilder>("empty") };
    DisplacementBuilder<UtilityIntegral> disp(crosss);
}

TEST(constructor, interpolant)
{
    auto cross = std::make_shared<CrossSectionBuilder>("empty");
    cross->SetdEdx_function([](double energy) { return 1; });

    CrossSectionList crosss{ cross };

    DisplacementBuilder<UtilityInterpolant> disp(crosss);
}

TEST(constructor, error)
{
    CrossSectionList crosss{};
    ASSERT_THROW(DisplacementBuilder<UtilityIntegral> disp(crosss),
        std::invalid_argument);
}

constexpr const std::array<double, 7> energies
    = { 1e2, 1e4, 1e6, 1e8, 1e10, 1e12, 1e14 };

TEST(utility_integral, function_to_integral)
{
    auto const_loss = std::make_shared<CrossSectionBuilder>("constant");
    const_loss->SetdEdx_function([](double energy) { return a; });
    auto lin_loss = std::make_shared<CrossSectionBuilder>("lin");
    lin_loss->SetdEdx_function([](double energy) { return b * energy; });

    CrossSectionList crosss{ const_loss, lin_loss };
    auto disp = std::make_shared<DisplacementBuilder<UtilityIntegral>>(crosss);

    for (const auto e : energies) {
        auto calc_result = disp->FunctionToIntegral(e);
        ASSERT_EQ(calc_result, -1 / (a + b * e));
    }
}

TEST(utility_integral, solve_track_integral)
{
    auto const_loss = std::make_shared<CrossSectionBuilder>("constant");
    const_loss->SetdEdx_function([](double energy) { return a; });
    auto lin_loss = std::make_shared<CrossSectionBuilder>("lin");
    lin_loss->SetdEdx_function([](double energy) { return b * energy; });

    CrossSectionList crosss{ const_loss, lin_loss };
    auto disp = std::make_shared<DisplacementBuilder<UtilityIntegral>>(crosss);

    std::function<double(double, double)> displacement_analytic =
        [](double energy_initial, double energy_final) {
            return (std::log((a + b * energy_initial) / (a + b * energy_final)))
                / b;
        };

    for (size_t i = 0; i < energies.size() - 1; ++i) {
        auto disp_calc
            = disp->SolveTrackIntegral(energies[i + 1], energies[i], 0);
        auto disp_analytic
            = displacement_analytic(energies[i + 1], energies[i]);

        EXPECT_NEAR(disp_calc, disp_analytic, PARTICLE_POSITION_RESOLUTION);
    }
}

TEST(utility_integral, upper_limit_track_integral)
{
    // dEdx = 1 + E
    auto const_loss = std::make_shared<CrossSectionBuilder>("constant");
    const_loss->SetdEdx_function([](double energy) { return a; });
    auto lin_loss = std::make_shared<CrossSectionBuilder>("lin");
    lin_loss->SetdEdx_function([](double energy) { return b * energy; });

    CrossSectionList crosss{ const_loss, lin_loss };
    auto disp = std::make_shared<DisplacementBuilder<UtilityIntegral>>(crosss);

    std::function<double(double, double)> upper_limit_analytic
        = [](double energy_initial, double disp) {
              return ((a + b * energy_initial) * std::exp(-b * disp) - a) / b;
          };

    for (const auto& e : energies) {
        auto energy_calc = disp->UpperLimitTrackIntegral(e, 100);
        auto energy_analytic = upper_limit_analytic(e, 100);

        EXPECT_NEAR(energy_calc, energy_analytic, energy_calc * PARTICLE_POSITION_RESOLUTION);
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
