#pragma once
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/propagation_utility/ContRand.h"
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
namespace PROPOSAL {
template <class T> class ContRandBuilder : public ContRand {
    T cont_rand_integral;

    void build_tables() {};

public:
    template <typename Cross>
    ContRandBuilder(std::shared_ptr<Displacement> disp, Cross&& cross)
        : ContRand(disp, std::forward<Cross>(cross))
        , cont_rand_integral([this](double E) { return FunctionToIntegral(E); },
              disp->GetLowerLim(), this->GetHash())
    {
    }

    inline double Variance(double initial_energy, double final_energy) final
    {
        assert(initial_energy >= final_energy);
        return cont_rand_integral.Calculate(initial_energy, final_energy);
    }

    inline double EnergyRandomize(
        double initial_energy, double final_energy, double rnd) final
    {
        auto std = std::sqrt(Variance(initial_energy, final_energy));
        return SampleFromGaussian(
            final_energy, std, rnd, disp->GetLowerLim(), initial_energy);
    }
};

template <typename T>
auto make_contrand(
    std::shared_ptr<Displacement> disp, T&& cross, bool interpolate = true)
{
    auto cont_rand = std::unique_ptr<ContRand>();
    if (interpolate)
        cont_rand = std::make_unique<ContRandBuilder<UtilityInterpolant>>(
            disp, std::forward<T>(cross));
    else
        cont_rand = std::make_unique<ContRandBuilder<UtilityIntegral>>(
            disp, std::forward<T>(cross));
    return cont_rand;
}

template <typename T> auto make_contrand(T&& cross, bool interpolate = true)
{
    auto disp = std::shared_ptr<Displacement>(make_displacement(cross, false));
    return make_contrand(disp, cross, interpolate);
}
} // namespace PROPOSAL
