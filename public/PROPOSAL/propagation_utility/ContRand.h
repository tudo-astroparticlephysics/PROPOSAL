#pragma once
#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/particle/ParticleDef.h"

namespace PROPOSAL {

class ContRand {
    public:
        ContRand(CrossSectionList cross, const ParticleDef &def) : cross(cross), mass(def.mass) {}
        virtual double EnergyRandomize(double initial_energy, double final_energy, double rnd) = 0;

    protected:
        CrossSectionList cross;
        double mass;
        std::string name = "contrand";
    };

    template<class T>
class ContRandBuilder : public ContRand {
    public:
        ContRandBuilder<T>(CrossSectionList cross, const ParticleDef &def)
                : ContRand(cross, def),
                  displacement(cross),
                  integral(std::bind(&ContRandBuilder::ContRandIntegrand, this, std::placeholders::_1)) {
            if(typeid(T) == typeid(UtilityInterpolant)){
                size_t hash_digest = 0;
                for (const auto& crosssection: cross){
                    hash_combine(hash_digest, crosssection->GetHash());
                }
                integral.BuildTables(name, hash_digest);
            }
        }

        double ContRandIntegrand (double energy) {
            assert(energy >= 0);
            double sum = 0.0;
            for (const auto &crosssections : cross)
                sum += crosssections->CalculatedE2dx(energy);

            return displacement.FunctionToIntegral(energy) * sum;
        }

        double EnergyRandomize(double initial_energy, double final_energy, double rnd) override {
            assert(initial_energy >= final_energy);
            double variance = integral.Calculate(initial_energy, final_energy, 0.0);
            return SampleFromGaussian(final_energy, std::sqrt(variance), rnd, mass, initial_energy);
        }

    private:
        T integral;
        DisplacementBuilder <UtilityIntegral> displacement;

    };

}