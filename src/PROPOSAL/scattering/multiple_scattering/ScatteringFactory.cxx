#include "PROPOSAL/scattering/multiple_scattering/ScatteringFactory.h"

namespace PROPOSAL {
namespace detail {
    std::unique_ptr<multiple_scattering::Parametrization>
    make_multiple_scattering(
        ScatteringType type, ParticleDef const& p_def, Medium const& medium)
    {
        switch (type) {
        case ScatteringType::Highland:
            return unique_ptr<multiple_scattering::Parametrization>(
                new multiple_scattering::Highland(p_def, medium));
        case ScatteringType::Moliere:
            return unique_ptr<multiple_scattering::Parametrization>(
                new multiple_scattering::Moliere(p_def, medium));
        case ScatteringType::NoScattering:
            return nullptr;
        default:
            throw std::out_of_range("This constructor is not provided.");
        }
    }
} // namespace detail
} // namespace PROPOSAL
