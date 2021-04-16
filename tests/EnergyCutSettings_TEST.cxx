#include "gtest/gtest.h"

#include "PROPOSAL/EnergyCutSettings.h"
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <vector>
#include <numeric>

using nlohmann::json;
using std::string;
using std::vector;
using namespace PROPOSAL;

TEST(constructor, parameters)
{
    const auto ecut = 123;
    const auto vcut = 0.123;
    const auto cont_rand = true;

    EnergyCutSettings cut(ecut, vcut, cont_rand);

    EXPECT_EQ(cut.GetEcut(), ecut) << "Energy cut initalization failed.";
    EXPECT_EQ(cut.GetVcut(), vcut) << "Relative cut initalization failed.";
    EXPECT_EQ(cut.GetContRand(), cont_rand)
        << "Continuous randomization initalization.";
}

TEST(constructor, json)
{
    auto config = json::parse("{ \"e_cut\": 123, \"v_cut\": 0.123, \"cont_rand\": true }");
    EnergyCutSettings cut(config);

    EXPECT_EQ(cut.GetEcut(), 123) << "Energy cut initalization failed.";
    EXPECT_EQ(cut.GetVcut(), 0.123) << "Relative cut initalization failed.";
    EXPECT_EQ(cut.GetContRand(), true)
        << "Continuous randomization initalization.";
}

TEST(constructor, copy)
{
    const auto ecut = 123;
    const auto vcut = 0.123;
    const auto cont_rand = true;

    EnergyCutSettings cut1(ecut, vcut, cont_rand);
    EnergyCutSettings cut2(cut1);

    EXPECT_EQ(cut2.GetEcut(), ecut) << "Energy cut initalization failed.";
    EXPECT_EQ(cut2.GetVcut(), vcut) << "Relative cut initalization failed.";
    EXPECT_EQ(cut2.GetContRand(), cont_rand) << "Continuous randomization initalization.";
}

TEST(constructor, exceptions)
{
    ASSERT_THROW(EnergyCutSettings(0, -1.0, true), std::invalid_argument);
    ASSERT_THROW(EnergyCutSettings(0, 1.1, true), std::invalid_argument);
    ASSERT_THROW(EnergyCutSettings(-500, 0.05, true), std::invalid_argument);
}

TEST(operator, comparison)
{
    EnergyCutSettings cut1(123, 0.123, true);
    EnergyCutSettings cut2(cut1);

    ASSERT_EQ(cut1, cut2);
}

TEST(cut, absolut)
{
    vector<double> energies{ 1, 1e3, 1e6, 1e9, 1e12 };
    auto vcut = 1.0;
    auto ecut = 500;

    EnergyCutSettings cut(ecut, vcut, true);
    for (const auto& e : energies)
        EXPECT_EQ((ecut / e < vcut) ? ecut / e : vcut, cut.GetCut(e))
            << "wrong energy cut behaviour.";
}

TEST(cut, relative)
{
    vector<double> energies{ 1, 1e3, 1e6, 1e9, 1e12 };
    auto vcut = 0.05;
    auto ecut = std::numeric_limits<double>::infinity();

    EnergyCutSettings cut(ecut, vcut, true);
    for (const auto& e : energies)
        EXPECT_EQ((e * vcut < ecut) ? vcut : ecut / e, cut.GetCut(e))
            << "wrong relative cut behaviour.";
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
