#include "PROPOSAL/scattering/multiple_scattering/HighlandIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"

using namespace PROPOSAL;
using namespace multiple_scattering;

double HighlandIntegral::CalculateTheta0(double grammage, double ei, double ef) {
    auto integral_result = highland_integral->Calculate(ei, ef);
    assert(integral_result >= 0);
    auto aux = 13.6
               * std::sqrt(std::max(integral_result, 0.) / radiation_length)
               * std::abs(charge);
    aux *= std::max(
            1. + 0.038 * std::log(grammage / radiation_length), 0.0);
    return std::min(aux, 1.0);
}

double HighlandIntegral::Integral(Displacement& disp, double energy)
{
    auto square_momentum = (energy - mass) * (energy + mass);
    auto aux = energy / square_momentum;
    return disp.FunctionToIntegral(energy) * aux * aux;
}

HighlandIntegral::HighlandIntegral(const ParticleDef& p, Medium const& m,
                 std::shared_ptr<Displacement> disp, std::false_type)
        : Highland(p, m)
        , highland_integral(std::make_unique<UtilityIntegral>(
                [this, disp](double E) { return Integral(*disp, E); },
                disp->GetLowerLim(), disp->GetHash())) {};

HighlandIntegral::HighlandIntegral(const ParticleDef& p, Medium const& m,
                 std::shared_ptr<Displacement> disp, std::true_type)
        : Highland(p, m)
        , highland_integral(std::make_unique<UtilityInterpolant>(
                [this, disp](double E) { return Integral(*disp, E); },
                disp->GetLowerLim(), disp->GetHash())) {
        highland_integral->BuildTables("scattering_", 500, true);
};

namespace PROPOSAL {
    std::unique_ptr<multiple_scattering::Parametrization> make_highland_integral(
            ParticleDef const &p, Medium const &m,
            std::shared_ptr<Displacement> disp, bool interpol) {
        auto scatter = std::unique_ptr<multiple_scattering::Parametrization>();
        if (interpol)
            scatter = std::make_unique<multiple_scattering::HighlandIntegral>(
                    p, m, disp, std::true_type{});
        else
            scatter = std::make_unique<multiple_scattering::HighlandIntegral>(
                    p, m, disp, std::false_type{});
        return scatter;
    }

    std::unique_ptr<multiple_scattering::Parametrization> make_highland_integral(
            ParticleDef const &p, Medium const &m,
            std::vector<std::shared_ptr<CrossSectionBase>> const &c,
            bool interpol) {
        auto disp = std::shared_ptr<Displacement>(make_displacement(c, false));
        return make_highland_integral(p, m, disp, interpol);
    }
}