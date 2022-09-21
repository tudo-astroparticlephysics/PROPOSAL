#include "PROPOSAL/propagation_utility/ContRand.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/propagation_utility/Displacement.h"

using namespace PROPOSAL;

Interpolant1DBuilder::Definition ContRand::interpol_def = { 200 };

ContRand::ContRand(std::shared_ptr<Displacement> _disp, std::vector<std::shared_ptr<CrossSectionBase>> const& cross)
    : disp(_disp)
    , hash(0)
    , cross_list(cross)
{
    // Check if there is at least one crosssection where dE2dx tables have been built.
    // Otherwise, it doesn't make sense to build a dE2dx object, so we throw an exception.
    bool has_dE2dx = false; // is there at least one crosssection with cont_rand enabled?
    for (auto& c : cross) {
        if (c->GetEnergyCutSettings()) // avoid accessing nullptr
            has_dE2dx = has_dE2dx || c->GetEnergyCutSettings()->GetContRand();
    }
    if (has_dE2dx == false)
        throw std::invalid_argument("Can not build ContRand because no dE2dx tables have been built. You need to use "
                                    "continuous_randomization=True when creating the EnergyCutSettings object to be "
                                    "able to use continuous randomization.");
    CalculateHash();
}

double ContRand::FunctionToIntegral(double energy)
{
    assert(energy >= 0);
    double sum = 0.0;
    for (auto& crosssections : cross_list)
        sum += crosssections->CalculatedE2dx(energy);
    return disp->FunctionToIntegral(energy) * sum;
}

void ContRand::CalculateHash() noexcept {
    hash_combine(hash, disp->GetHash());
}