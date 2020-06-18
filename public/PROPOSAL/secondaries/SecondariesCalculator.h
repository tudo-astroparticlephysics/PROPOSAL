#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/secondaries/Parametrization.h"

#include <memory>
#include <unordered_map>
#include <vector>

using std::unordered_map;
using std::vector;

namespace PROPOSAL {
class SecondariesCalculator {
    using param_ptr = unique_ptr<secondaries::Parametrization>;
    unordered_map<int, param_ptr> secondary_generator;

public:
    SecondariesCalculator() = default;

    template <typename Param> void addInteraction(Param&&);

    vector<Loss::secondary_t> CalculateSecondaries(double primary_energy,
        Loss::secondary_t loss, const Component& comp, vector<double> rnd);

    size_t requiredRandomNumbers(int);
};
}

namespace PROPOSAL {
template <typename Param>
void SecondariesCalculator::addInteraction(Param&& param)
{
    secondary_generator[static_cast<int>(param.type)]
        = param_ptr(make_unique<typename std::decay<Param>::type>(std::forward<Param>(param)));
}
}
