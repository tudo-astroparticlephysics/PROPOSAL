#pragma once

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/propagation_utility/Displacement.h"

namespace PROPOSAL {

class Time {
public:
    Time();
    virtual double TimeElapsed(
        double initial_energy, double final_energy, double time)
        = 0;
    virtual double TimeElapsed(double distance) = 0;

protected:
    std::string name = "time";
    double lower_lim;
};

extern Interpolant1DBuilder::Definition time_interpol_def;

template <class T> class ExactTimeBuilder : public Time {
public:
    ExactTimeBuilder<T>(CrossSectionList cross, const ParticleDef& def)
        : mass(cross.front()->GetParametrization().GetParticleMass())
        , displacement(cross)
        , integral(std::bind(
              &ExactTimeBuilder::TimeIntegrand, this, std::placeholders::_1), mass)
    {
        if (typeid(T) == typeid(UtilityInterpolant)) {
            size_t hash_digest = 0;
            for (const auto& c : cross)
                hash_combine(hash_digest, c->GetParametrization().GetHash(),
                    c->GetParametrization().GetMultiplier());
            integral.BuildTables(name, hash_digest, time_interpol_def);
        }
    }

    double TimeIntegrand(double energy)
    {
        assert(energy > mass);

        auto square_momentum = (energy - mass) * (energy + mass);
        auto particle_momentum = std::sqrt(square_momentum);

        return displacement.FunctionToIntegral(energy) * energy
            / (particle_momentum * SPEED);
    }

    double TimeElapsed(
        double initial_energy, double final_energy, double time) override
    {
        return integral.Calculate(initial_energy, final_energy, time);
    }

    double TimeElapsedUpperLimit(double initial_energy, double time)
    {
        return integral.GetUpperLimit(initial_energy, time);
    }

    double TimeElapsed(double distance) override
    {
        throw std::logic_error(
            "Exact elapsed time can only be calculated using two energies");
    }


private:
    double mass;
    DisplacementBuilder<UtilityIntegral> displacement;
    T integral;
};

class ApproximateTimeBuilder : public Time {
public:
    ApproximateTimeBuilder() {}

    double TimeElapsed(
        double initial_energy, double final_energy, double time) override
    {
        (void)initial_energy;
        (void)final_energy;
        (void)time;
        throw std::logic_error("Appoximated elapsed time can only be "
                               "calculated using a given distance");
    }

    double TimeElapsed(double distance) override { return distance / SPEED; }
};
template <class T>
Interpolant1DBuilder::Definition ExactTimeBuilder<T>::time_interpol_def;
}
