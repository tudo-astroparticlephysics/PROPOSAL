#pragma once

#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include <memory>
#include <spdlog/fwd.h>

namespace PROPOSAL {
class EnergyCutSettings;
enum class InteractionType;
} // namespace PROPOSAL

/* namespace PROPOSAL { */
/* namespace detail { */
/*     size_t dEdx_Hash(size_t, crosssection::Parametrization<Medium> const&, */
/*         ParticleDef const&, Medium const&, EnergyCutSettings const&); */

/*     size_t dEdx_Hash(size_t, crosssection::Parametrization<Component> const&,
 */
/*         ParticleDef const&, Component const&, EnergyCutSettings const&); */
/* } // namespace detail */
/* } // namespace PROPOSAL */

namespace PROPOSAL {
class CrossSectionDEDX {
protected:
    size_t hash;
    std::shared_ptr<spdlog::logger> logger;

public:
    CrossSectionDEDX(size_t hash);

    virtual ~CrossSectionDEDX() = default;

    virtual double Calculate(double energy) const = 0;

    size_t GetHash() const noexcept { return hash; }
};
} // namespace PROPOSAL
