/* #pragma once */
/* #include "PROPOSAL/math/MathMethods.h" */
/* #include "PROPOSAL/particle/ParticleDef.h" */
/* #include "PROPOSAL/Constants.h" */
/* #include "PROPOSAL/propagation_utility/Displacement.h" */
/* #include "PROPOSAL/math/InterpolantBuilder.h" */
/* #include "PROPOSAL/crossection/CrossSection.h" */

/* namespace PROPOSAL { */

/* class Decay { */
/* public: */
/*     Decay(const CrossSectionList& cross, const ParticleDef&); */
/*     Decay(const CrossSectionList& cross, double lifetime, double mass); */
/*     virtual ~Decay() {}; */
/*     virtual double EnergyDecay(double initial_energy, double rnd) = 0; */

/* protected: */
/*     CrossSectionList cross; */
/*     double mass; */
/*     double lifetime; */
/*     double lower_lim; */
/* }; */

/* extern Interpolant1DBuilder::Definition decay_interpol_def; */

/* template <class T> class DecayBuilder : public Decay { */
/* public: */
/*     DecayBuilder<T>(CrossSectionList cross, const ParticleDef& p_def) : DecayBuilder<T>(cross, p_def.lifetime, p_def.mass){}; */

/*     DecayBuilder<T>(CrossSectionList cross, double lifetime, double mass) */
/*         : Decay(cross, lifetime, mass) */
/*         , displacement(cross) */
/*         , integral(std::bind( */
/*               &DecayBuilder::DecayIntegrand, this, std::placeholders::_1), lower_lim) */
/*     { */
/*         if (typeid(T) == typeid(UtilityInterpolant)) { */
/*             size_t hash_digest = 0; */
/*             for (const auto& c: cross) */
/*                 hash_combine(hash_digest, c->GetHash()); */
/*             decay_interpol_def.function1d = [this](double energy) { */
/*                 return reinterpret_cast<UtilityIntegral*>(&integral)->Calculate( */
/*                         lower_lim, energy); */
/*             }; */
/*             integral.BuildTables("decay", hash_digest, decay_interpol_def); */
/*         } */
/*     } */

/*     double DecayIntegrand(double energy) */
/*     { */
/*         assert(lifetime < 0); */
/*         assert(energy <= mass); */

/*         double square_momentum = (energy - mass) * (energy + mass); */
/*         double aux = SPEED * std::sqrt(square_momentum) / mass; */
/*         return displacement.FunctionToIntegral(energy) / aux; */
/*     } */

/*     double EnergyDecay(double initial_energy, double rnd) override */
/*     { */
/*         auto rndd = -std::log(rnd); */
/*         auto rnddMin = integral.Calculate(initial_energy, lower_lim) / lifetime; */

/*         if (rndd >= rnddMin) */
/*             return lower_lim; */

/*         return integral.GetUpperLimit(initial_energy, rndd) / lifetime; */
/*     } */


/* private: */
/*     T integral; */
/*     DisplacementBuilder<UtilityIntegral> displacement; */
/* }; */
/* } */
