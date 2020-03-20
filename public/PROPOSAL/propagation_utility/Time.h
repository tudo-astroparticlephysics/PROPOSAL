#pragma once
#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/propagation_utility/Displacement.h"
#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/ParticleDef.h"

namespace PROPOSAL {

class Time {
    public:
        Time() {}
        virtual double ElapsedTime(double initial_energy, double final_energy, double time) = 0;
        virtual double ElapsedTime(double distance) = 0;

    protected:
        std::string name = "time";
    };

    template<class T>
class ExactTimeBuilder : public Time {
    public:
        ExactTimeBuilder<T>(CrossSectionList cross, const ParticleDef &def)
            : mass(def.mass),
              displacement(cross),
              integral (std::bind( &ExactTimeBuilder::TimeIntegrand, this, std::placeholders::_1)){
            if(typeid(T) == typeid(UtilityInterpolant)){
                size_t hash_digest = 0;
                for (const auto& crosssection: cross){
                    hash_combine(hash_digest, crosssection->GetParametrization().GetHash(),
                                 crosssection->GetParametrization().GetMultiplier());
                }
                integral.BuildTables(name, hash_digest);
            }
        }

        double TimeIntegrand(double energy) {
            double square_momentum{(energy - mass) * (energy + mass)};
            double particle_momentum{std::sqrt(std::max(square_momentum, 0.0))};
            return displacement.FunctionToIntegral(energy) * energy / (particle_momentum * SPEED);
        }

        double ElapsedTime(double initial_energy, double final_energy, double time) override {
            return integral.Calculate(initial_energy, final_energy, time);
        }

        double ElapsedTimeUpperLimit(double initial_energy, double time) {
            return integral.GetUpperLimit(initial_energy, time);
        }

        double ElapsedTime(double distance) override {
            throw std::logic_error("Exact elapsed time can only be calculated using two energies");
        }

    private:
        CrossSectionList cross;
        T integral;
        double mass;
        DisplacementBuilder<UtilityIntegral> displacement;
    };

class ApproximateTimeBuilder : public Time {
    public:
        ApproximateTimeBuilder() {}

        double ElapsedTime(double initial_energy, double final_energy, double time) override {
            (void)initial_energy;
            (void)final_energy;
            (void)time;
            throw std::logic_error("Appoximated elapsed time can only be calculated using a given distance");
        }

        double ElapsedTime(double distance) override { return distance / SPEED; }
    };

}
