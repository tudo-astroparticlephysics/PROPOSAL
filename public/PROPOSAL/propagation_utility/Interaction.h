#pragma once
#include <tuple>
#include "PROPOSAL/propagation_utility/Displacement.h"

using std::tuple;

using comp_rates = unordered_map<size_t, double>;

namespace PROPOSAL {

class Interaction {
    CrossSectionList cross;
    DisplacementBuilder<UtilityIntegral> displacement;

protected:
    double lower_lim;

    double FunctionToIntegral(double);

public:
    Interaction(CrossSectionList cross);
    virtual ~Interaction() = default;


    unordered_map<InteractionType, comp_rates> CalculateRates(double);

    /* enum { RELATIV_LOSS, INTERACTION_TYPE }; */
    tuple<double, InteractionType> SampleStochasticLoss(double, double);

    virtual double SolveEnergyIntegral(double, double) = 0;
    virtual double UpperLimitEnergyIntegral(double, double) = 0;
};

extern Interpolant1DBuilder::Definition interaction_interpol_def;

template <class T> class InteractionBuilder : public Interaction {
    T integral;

public:
    InteractionBuilder<T>(CrossSectionList cross)
        : Interaction(cross)
        , integral(
              std::bind(&Interaction::FunctionToIntegral, this, _1), lower_lim)
    {
        if (typeid(T) == typeid(UtilityInterpolant)) {
            size_t hash_digest = 0;
            for (const auto& c : cross)
                hash_combine(hash_digest, c->GetHash());
            interaction_interpol_def.function1d = [this](double energy) {
                return reinterpret_cast<UtilityIntegral*>(&integral)->Calculate(
                    energy, 1e14);
            };
            integral.BuildTables(
                "interaction", hash_digest, interaction_interpol_def);
        }
    }

    double SolveEnergyIntegral(
        double initial_energy, double final_energy) override
    {
        return integral.Calculate(initial_energy, final_energy);
    }

    double UpperLimitEnergyIntegral(double initial_energy, double rnd) override
    {
        return integral.GetUpperLimit(initial_energy, -std::log(rnd));
    }
};
}
