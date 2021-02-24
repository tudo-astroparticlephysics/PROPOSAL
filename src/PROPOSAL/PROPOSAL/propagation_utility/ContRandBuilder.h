#pragma once
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/propagation_utility/ContRand.h"
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
namespace PROPOSAL {

template <class T> class ContRandBuilder : public ContRand {
    T cont_rand_integral;

    void build_tables();

public:
    ContRandBuilder(std::shared_ptr<Displacement> disp,
        std::vector<std::shared_ptr<CrossSectionBase>> const& cross)
        : ContRand(disp, cross)
        , cont_rand_integral([this](double E) { return FunctionToIntegral(E); },
              disp->GetLowerLim(), this->GetHash())
    {
        build_tables();
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
        return SampleFromGaussian(final_energy, std, rnd,
                                  std::min(disp->GetLowerLim(), final_energy),
                                  initial_energy);
    }
};

std::unique_ptr<ContRand> make_contrand(std::shared_ptr<Displacement>,
    std::vector<std::shared_ptr<CrossSectionBase>> const&,
    bool interpolate = true);

std::unique_ptr<ContRand> make_contrand(
    std::vector<std::shared_ptr<CrossSectionBase>> const&,
    bool interpolate = true);
} // namespace PROPOSAL
