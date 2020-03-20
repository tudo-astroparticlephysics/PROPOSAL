#pragma once
#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/particle/ParticleDef.h"


namespace PROPOSAL {

class Decay {
    public:
        Decay(CrossSectionList cross, const ParticleDef &def) : cross(cross), mass(def.mass), lifetime(def.lifetime) {};
        virtual double EnergyDecay(double initial_energy, double rnd) = 0;

    protected:
        CrossSectionList cross;
        double mass;
        double lifetime;
        std::string name = "decay";

    };

    template<class T>
class DecayBuilder : public Decay {
    public:
        DecayBuilder<T>(CrossSectionList cross, const ParticleDef &def)
            : Decay(cross, def),
              displacement(cross),
              integral(std::bind(&DecayBuilder::DecayIntegrand, this, std::placeholders::_1)){
            if(typeid(T) == typeid(UtilityInterpolant)){
                size_t hash_digest = 0;
                for (const auto& crosssection: cross){
                    hash_combine(hash_digest, crosssection->GetParametrization().GetHash(),
                                 crosssection->GetParametrization().GetMultiplier());
                }
                integral.BuildTables(name, hash_digest);
            }
        }

        double DecayIntegrand (double energy) {
            if (lifetime < 0)
                return 0;
            double square_momentum = (energy - mass) * (energy + mass);
            double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));
            double aux = 1.0 / std::max(SPEED * particle_momentum / mass, PARTICLE_POSITION_RESOLUTION);
            return displacement.FunctionToIntegral(energy) * aux;
        }

        double EnergyDecay(double initial_energy, double rnd) override {
            auto rndd = -std::log(rnd);
            auto rnddMin = 0;
            rnddMin = integral.Calculate(initial_energy, mass, rndd) / lifetime;

            if (rndd >= rnddMin || rnddMin <= 0)
                return mass;

            return integral.GetUpperLimit(initial_energy, rndd) / lifetime;
        }

    private:
        T integral;
        DisplacementBuilder <UtilityIntegral> displacement;
    };

}