#pragma once

#include <memory>
#include <spdlog/fwd.h>

namespace PROPOSAL {
namespace crosssection {
    template <typename T> class Parametrization;
} // namespace crosssection
class ParticleDef;
class Medium;
class Component;
class EnergyCutSettings;
enum class InteractionType;
} // namespace PROPOSAL

namespace PROPOSAL {
class CrossSectionDEDX {
protected:
    size_t hash;
    std::shared_ptr<spdlog::logger> logger;

public:
    CrossSectionDEDX(size_t hash, std::string param_name);

    virtual ~CrossSectionDEDX() = default;

    virtual double Calculate(double energy) const = 0;

    size_t GetHash() const noexcept { return hash; }
};
} // namespace PROPOSAL

namespace PROPOSAL {
namespace detail {
    size_t dEdx_Hash(InteractionType,
        crosssection::Parametrization<Medium> const&, ParticleDef const&,
        Medium const&, EnergyCutSettings const&);

    size_t dEdx_Hash(InteractionType,
        crosssection::Parametrization<Component> const&, ParticleDef const&,
        Component const&, EnergyCutSettings const&);
} // namespace detail
} // namespace PROPOSAL
