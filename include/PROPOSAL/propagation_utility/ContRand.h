#pragma once

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/CrossSectionVector.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/propagation_utility/Displacement.h"

namespace PROPOSAL {

struct ContRand {
    using crossbase_list_t = std::vector<std::shared_ptr<CrossSectionBase>>;

    std::shared_ptr<Displacement> disp;
    size_t hash;

    crossbase_list_t cross_list;

    template <typename Cross>
    ContRand(std::shared_ptr<Displacement> _disp, Cross&& cross)
        : disp(_disp)
        , hash(0)
        , cross_list(std::begin(cross), std::end(cross))
    {
        hash_combine(hash, _disp->GetHash());
    }

    virtual ~ContRand() = default;

    static Interpolant1DBuilder::Definition interpol_def;

    double FunctionToIntegral(double);

    virtual double Variance(double, double) = 0;
    virtual double EnergyRandomize(double, double, double) = 0;

    auto GetHash() const noexcept { return hash; }
};

} // namespace PROPOSAL
