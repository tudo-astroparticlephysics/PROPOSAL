#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/propagation_utility/ContRand.h"
#include "PROPOSAL/math/MathMethods.h"

namespace PROPOSAL {
template <class T> class ContRandBuilder : public ContRand {
    T cont_rand_integral;

    T BuildContRandIntegral(crossbase_list_t const& cross)
    {
        auto disp = std::shared_ptr<Displacement>(make_displacement(cross, false));
        auto cont_rand_func = [this, disp](double energy) mutable {
            return FunctionToIntegral(*disp, energy);
        };
        T cont_rand_integral(cont_rand_func, lower_lim);
        if (typeid(T) == typeid(UtilityInterpolant)) {
            auto hash_digest = (size_t)0;
            hash_combine(hash_digest, CrossSectionVector::GetHash(cross), interpol_def->GetHash());

            if(not interpol_def)
                interpol_def = std::make_unique<Interpolant1DBuilder::Definition>();

            cont_rand_integral.BuildTables("contrand", hash_digest, *interpol_def);
        };
        return cont_rand_integral;
    }

    public:
    template <typename Cross>
    ContRandBuilder(Cross&& cross)
        : ContRand(std::forward<Cross>(cross))
        , cont_rand_integral(BuildContRandIntegral(cross_list))
    {
    }

    inline double Variance(double initial_energy, double final_energy) final
    {
        assert(initial_energy >= final_energy);
        return cont_rand_integral.Calculate(initial_energy, final_energy);
    }

    inline double EnergyRandomize(
            double initial_energy, double final_energy, double rnd) final
    {
        auto std = std::sqrt(Variance(initial_energy, final_energy));
        return SampleFromGaussian(final_energy, std, rnd, lower_lim, initial_energy);
    }
};

template <typename T>
std::unique_ptr<ContRand> make_contrand(T&& cross, bool interpolate = true)
{
    if (interpolate)
        return PROPOSAL::make_unique<ContRandBuilder<UtilityInterpolant>>(
                std::forward<T>(cross));
    return PROPOSAL::make_unique<ContRandBuilder<UtilityIntegral>>(
            std::forward<T>(cross));
}
} // namespace PROPOSAL
