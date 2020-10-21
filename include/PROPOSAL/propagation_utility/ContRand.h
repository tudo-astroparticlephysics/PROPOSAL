#pragma once

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/CrossSectionVector.h"
#include "PROPOSAL/propagation_utility/Displacement.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

namespace PROPOSAL {

struct ContRand {
    using crossbase_list_t = std::vector<std::shared_ptr<CrossSectionBase>>;

    crossbase_list_t cross_list;
    double lower_lim;

    template <typename Cross>
    ContRand(Cross&& cross)
            : cross_list(std::begin(cross), std::end(cross))
            , lower_lim(CrossSectionVector::GetLowerLim(cross))
    {}
    virtual ~ContRand() = default;

    static std::unique_ptr<Interpolant1DBuilder::Definition> interpol_def;

    double FunctionToIntegral(Displacement&, double);

    virtual double Variance(double, double) = 0;
    virtual double EnergyRandomize(double, double, double) = 0;
};

} // namespace PROPOSAL
