#pragma once
#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/Logging.h"

namespace PROPOSAL {

class Interaction {
    public:
        Interaction(CrossSectionList cross, const ParticleDef &def) : cross(cross), mass(def.mass) {}
        virtual double EnergyInteraction(double initial_energy, double rnd) = 0;
        std::shared_ptr<CrossSection> TypeInteraction(double energy, const std::array<double, 2>& rnd){
            std::vector<double> rates;
            for (const auto& crosssection : cross)
                rates.push_back(crosssection->CalculatedNdx(energy, rnd[1]));

            double total_rate{ std::accumulate(rates.begin(), rates.end(), 0.0) };
            log_debug("Total rate = %f, total rate weighted = %f", total_rate, total_rate * rnd[0]);

            double rates_sum = 0;
            for (size_t i = 0; i < rates.size(); i++) {
                rates_sum += rates[i];
                if (rates_sum >= total_rate * rnd[0])
                    return cross.at(i);
            }

            throw std::logic_error("Something went wrong during the total rate calculation.");
        }

    protected:
        CrossSectionList cross;
        double mass;
        std::string name = "interaction";
    };

    template<class T>
class InteractionBuilder : public Interaction {
    public:
        InteractionBuilder<T>(CrossSectionList cross, const ParticleDef &def)
            : Interaction(cross, def),
              displacement(cross),
              integral(std::bind(&InteractionBuilder::InteractionIntegrand, this, std::placeholders::_1)) {
            if(typeid(T) == typeid(UtilityInterpolant)){
                size_t hash_digest = 0;
                for (const auto& crosssection: cross){
                    hash_combine(hash_digest, crosssection->GetParametrization().GetHash(),
                                 crosssection->GetParametrization().GetMultiplier());
                }
                integral.BuildTables(name, hash_digest, interaction_interpol_def);
            }
        }

        double InteractionIntegrand (double energy) {
            double total_rate = 0.0;
            for (const auto& crosssection : cross)
                total_rate += crosssection->CalculatedNdx(energy);

            return displacement.FunctionToIntegral(energy) * total_rate;
        }

        double EnergyInteraction (double initial_energy, double rnd) override {
            auto rndi = -std::log(rnd);
            auto rndiMin = 0.;

            rndiMin = integral.Calculate(initial_energy, mass, rndi);

            if (rndi >= rndiMin || rndiMin <= 0)
                return mass;

            return displacement.UpperLimitTrackIntegral(initial_energy, rndi);
        }

        static Interpolant1DBuilder::Definition interaction_interpol_def;
    private:
        T integral;
        DisplacementBuilder<UtilityIntegral> displacement;
    };

    template <class T>
    Interpolant1DBuilder::Definition InteractionBuilder<T>::interaction_interpol_def;

}
