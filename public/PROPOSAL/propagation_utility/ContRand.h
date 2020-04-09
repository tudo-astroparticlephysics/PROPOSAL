#pragma once

#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/propagation_utility/Displacement.h"

namespace PROPOSAL {

class ContRand {
public:
    ContRand(CrossSectionList cross);
    virtual double EnergyRandomize(
        double initial_energy, double final_energy, double rnd)
        = 0;

protected:
    CrossSectionList cross;
    double lower_lim;
};

extern Interpolant1DBuilder::Definition contrand_interpol_def;

template <class T> class ContRandBuilder : public ContRand {
public:
    ContRandBuilder<T>(CrossSectionList cross)
        : ContRand(cross)
        , displacement(cross)
        , integral(std::bind(&ContRandBuilder::ContRandIntegrand, this, std::placeholders::_1), lower_lim)
    {
        if (cross.size() < 1)
            throw std::invalid_argument("at least one crosssection is required.");

        if (typeid(T) == typeid(UtilityInterpolant)) {
            size_t hash_digest = 0;
            for (const auto& c : cross) {
                hash_combine(hash_digest, c->GetHash());
            }
            contrand_interpol_def.function1d = [this](double energy) {
                return reinterpret_cast<UtilityIntegral*>(&integral)->Calculate(
                        energy, lower_lim, 0);
            };
            integral.BuildTables("contrand", hash_digest, contrand_interpol_def);
        }
    }

    double ContRandIntegrand(double energy)
    {
        assert(energy >= 0);
        double sum = 0.0;
        for (const auto& crosssections : cross)
            sum += crosssections->CalculatedE2dx(energy);

        return displacement.FunctionToIntegral(energy) * sum;
    }

    double EnergyRandomize(
        double initial_energy, double final_energy, double rnd) override
    {
        assert(initial_energy >= final_energy);
        double variance = integral.Calculate(initial_energy, final_energy, 0.0);
        return SampleFromGaussian(
            final_energy, std::sqrt(variance), rnd, lower_lim, initial_energy);
    }

private:
    T integral;
    DisplacementBuilder<UtilityIntegral> displacement;
};

}
