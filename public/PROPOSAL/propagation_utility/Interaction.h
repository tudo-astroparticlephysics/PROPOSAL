#pragma once
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/propagation_utility/Displacement.h"

namespace PROPOSAL {

class Interaction {
protected:
    CrossSectionList cross;
    double mass;

public:
    Interaction(CrossSectionList cross);
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
              mass)
    {
        if (typeid(T) == typeid(UtilityInterpolant)) {
            size_t hash_digest = 0;
            for (const auto& c : cross)
                hash_combine(hash_digest, c->GetHash());
            integral.BuildTables("interaction", hash_digest, interaction_interpol_def);
        }
    }

    double InteractionIntegrand(double energy)
    {
        double total_rate = 0.0;
        for (const auto& crosssection : cross)
            total_rate += crosssection->CalculatedNdx(energy);

        return displacement.FunctionToIntegral(energy) * total_rate;
    }

    double EnergyInteraction(double initial_energy, double rnd) override
    {
        auto rndi = -std::log(rnd);
        auto rndiMin = 0.;

        rndiMin = integral.Calculate(initial_energy, mass, rndi);

        if (rndi >= rndiMin || rndiMin <= 0)
            return mass;

        return displacement.UpperLimitTrackIntegral(initial_energy, rndi);
    }


private:
    T integral;
    DisplacementBuilder<UtilityIntegral> displacement;
};

template <class T>
Interpolant1DBuilder::Definition
    InteractionBuilder<T>::interaction_interpol_def;
}
