#pragma once

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/propagation_utility/Displacement.h"

namespace PROPOSAL {

class Time {
public:
    Time();
    virtual double TimeElapsed(double initial_energy, double final_energy, double time) = 0;
    virtual double TimeElapsed(double distance) = 0;

protected:
    std::string name = "time";
};

extern Interpolant1DBuilder::Definition time_interpol_def;

template <class T> class ExactTimeBuilder : public Time {
public:
    ExactTimeBuilder<T>(CrossSectionList cross, const ParticleDef& p_def) : ExactTimeBuilder<T>(cross, p_def.mass){};

    ExactTimeBuilder<T>(CrossSectionList cross, double mass)
        : mass(mass)
        , displacement(cross)
        , lower_lim(InitializeLowerLim(cross))
        , integral(std::bind(
              &ExactTimeBuilder::TimeIntegrand, this, std::placeholders::_1), lower_lim)
    {
        if (cross.size() < 1)
            throw std::invalid_argument("at least one crosssection is required.");



        if (typeid(T) == typeid(UtilityInterpolant)) {
            size_t hash_digest = 0;
            for (const auto& c : cross)
                hash_combine(hash_digest, c->GetParametrization().GetHash(),
                    c->GetParametrization().GetMultiplier());
            time_interpol_def.function1d = [this](double energy) {
                return reinterpret_cast<UtilityIntegral*>(&integral)->Calculate(
                        lower_lim, energy, 0);
            };
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
        (void) distance;
        throw std::logic_error(
            "Exact elapsed time can only be calculated using two energies");
    }


protected:
    double InitializeLowerLim(CrossSectionList cross){
        double lower_lim_tmp = std::numeric_limits<double>::max();
        for (auto c : cross)
            lower_lim_tmp = std::min(lower_lim_tmp, c->GetParametrization().GetLowerEnergyLim());
        return lower_lim_tmp;
    }
    double lower_lim;
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
}
