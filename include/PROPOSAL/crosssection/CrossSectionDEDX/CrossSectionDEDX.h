#pragma once

#include <memory>
#include <spdlog/fwd.h>

namespace PROPOSAL {
namespace crosssection {
    template <typename Target> class Parametrization;
} // namespace crosssection

class ParticleDef;
class Medium;
class Component;
class EnergyCutSettings;
} // namespace PROPOSAL

namespace PROPOSAL {
class CrossSectionDEDX {
protected:
    size_t hash;
    std::shared_ptr<spdlog::logger> logger;

public:
    CrossSectionDEDX(crosssection::Parametrization<Medium> const&,
        ParticleDef const&, Medium const&, EnergyCutSettings const&);

    CrossSectionDEDX(crosssection::Parametrization<Component> const&,
        ParticleDef const&, Component const&, EnergyCutSettings const&);

    virtual ~CrossSectionDEDX() = default;

    virtual double Calculate(double energy) const = 0;

    size_t GetHash() const noexcept { return hash; }
};
} // namespace PROPOSAL
