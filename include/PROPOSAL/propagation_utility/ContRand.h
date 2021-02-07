#pragma once

#include "PROPOSAL/math/InterpolantBuilder.h"
#include <vector>

namespace PROPOSAL {
class Displacement;
struct CrossSectionBase;
}

namespace PROPOSAL {

struct ContRand {
    using crossbase_list_t = std::vector<std::shared_ptr<CrossSectionBase>>;

    std::shared_ptr<Displacement> disp;
    size_t hash;

    crossbase_list_t cross_list;

    ContRand(std::shared_ptr<Displacement> _disp,
        std::vector<std::shared_ptr<CrossSectionBase>> const& cross)
        : disp(_disp)
        , hash(0)
        , cross_list(cross)
    {
        CalculateHash();
    }

    virtual ~ContRand() = default;

    static Interpolant1DBuilder::Definition interpol_def;

    double FunctionToIntegral(double);

    virtual double Variance(double, double) = 0;
    virtual double EnergyRandomize(double, double, double) = 0;

    auto GetHash() const noexcept { return hash; }

private:
    void CalculateHash() noexcept;
};

} // namespace PROPOSAL
