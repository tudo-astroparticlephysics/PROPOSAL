#pragma once
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/propagation_utility/Displacement.h"

namespace PROPOSAL {

class Interaction {
protected:
    CrossSectionList cross;
    double lower_lim;

public:
    Interaction(CrossSectionList cross);
    virtual ~Interaction() = default;
    virtual double EnergyInteraction(double, double) = 0;
    std::shared_ptr<CrossSection> TypeInteraction(double, const std::array<double, 2>&);
};

extern Interpolant1DBuilder::Definition interaction_interpol_def;

template <class T> class InteractionBuilder : public Interaction {
public:
    InteractionBuilder<T>(CrossSectionList cross)
        : Interaction(cross)
        , displacement(cross)
        , integral(std::bind(&InteractionBuilder::InteractionIntegrand, this,
                       std::placeholders::_1),
              lower_lim)
    {
        if (typeid(T) == typeid(UtilityInterpolant)) {
            size_t hash_digest = 0;
            for (const auto& c : cross)
                hash_combine(hash_digest, c->GetHash());
            interaction_interpol_def.function1d = [this](double energy) {
                return reinterpret_cast<UtilityIntegral*>(&integral)->Calculate(
                        energy, lower_lim, 0);
            };
            integral.BuildTables("interaction", hash_digest, interaction_interpol_def);
        }
    }

    double InteractionIntegrand(double energy)
    {
        assert(energy >= lower_lim);
        double total_rate = 0.0;
        for (const auto& crosssection : cross)
            total_rate += crosssection->CalculatedNdx(energy);

        return displacement.FunctionToIntegral(energy) * total_rate;
    }

    double EnergyInteraction(double initial_energy, double rnd) override
    {
        assert(initial_energy >= lower_lim);
        auto rndi = -std::log(rnd);
        auto rndiMin = integral.Calculate(initial_energy, lower_lim, rndi);
        if (rndi >= rndiMin)
            return lower_lim;
        return integral.GetUpperLimit(initial_energy, rndi);
    }


private:
    DisplacementBuilder<UtilityIntegral> displacement;
    T integral;
};

}
