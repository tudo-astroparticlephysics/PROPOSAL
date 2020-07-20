#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/DefaultFactory.h"
#include "PROPOSAL/secondaries/Parametrization.h"

#include <memory>
#include <unordered_map>
#include <vector>

using std::vector;

namespace PROPOSAL {
class SecondariesCalculator {
    using param_ptr = unique_ptr<secondaries::Parametrization>;
    std::unordered_map<InteractionType, param_ptr, InteractionType_hash> secondary_generator;

public:
    SecondariesCalculator() = default;

    template <typename ParamList> SecondariesCalculator(ParamList&& param_list)
    {
        for (auto&& p : param_list)
            addInteraction(move(p));
    }

    template <typename TypeList>
    SecondariesCalculator(
        TypeList type_list, const ParticleDef& p, const Medium& m)
    {
       for (auto& t : type_list)
            addInteraction(secondaries::DefaultFactory::Create(t, p, m));
    }

    template <typename Param> void addInteraction(Param&& p)
    {
        secondary_generator[p->GetInteractionType()] = move(p);
    }

    size_t RequiredRandomNumbers(InteractionType) const noexcept;

    vector<Loss::secondary_t> CalculateSecondaries(
        double primary_energy,
        Loss::secondary_t loss,
        const Component& comp,
        vector<double> rnd);

    size_t requiredRandomNumbers(InteractionType);
};

template <typename TypeList>
std::unique_ptr<SecondariesCalculator> make_secondaries(TypeList&& list, const ParticleDef& p, const Medium& m)
{
    return PROPOSAL::make_unique<SecondariesCalculator>(std::forward<TypeList>(list), p, m);
}

}
