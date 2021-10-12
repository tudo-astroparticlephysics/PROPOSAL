#include "gtest/gtest.h"
#include <cmath>

#include "PROPOSAL/secondaries/parametrization/ionization/NaivIonization.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

TEST(NaivIonization, EnergyMomentumConservation)
{
    auto E_i = 1e2;

    auto param = secondaries::NaivIonization(EPlusDef(), StandardRock());
    auto loss = StochasticLoss((int)InteractionType::Ioniz, 1e1, Cartesian3D(0, 0, 0), Cartesian3D(0, 0, 1), 0, 0, E_i);

    auto rnd = std::vector<double>{0.1909};
    auto secs = param.CalculateSecondaries(loss, Components::StandardRock(), rnd);

    auto p_sum = Cartesian3D(0., 0., 0.);
    double E_sum = 0.;

    for (auto sec : secs) {
        p_sum = p_sum + sec.GetMomentum() * sec.direction;
        E_sum += sec.energy;
    }

    EXPECT_EQ(secs.size(), 2);

    // check momentum conservation
    EXPECT_NEAR(p_sum.GetX(), 0., COMPUTER_PRECISION);
    EXPECT_NEAR(p_sum.GetY(), 0., COMPUTER_PRECISION);
    EXPECT_NEAR(p_sum.GetZ(), std::sqrt((E_i + ME) * (E_i - ME)), COMPUTER_PRECISION);

    // check energy conservation
    EXPECT_NEAR(E_sum, E_i + ME, (E_i + ME) * COMPUTER_PRECISION);

}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
