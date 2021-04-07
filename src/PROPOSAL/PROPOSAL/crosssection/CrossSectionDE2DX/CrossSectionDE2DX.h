#pragma once

#include <memory>
#include <spdlog/fwd.h>

namespace PROPOSAL {
class EnergyCutSettings;
struct ParticleDef;
class Medium;
class Component;
namespace crosssection {
    template <typename T> class Parametrization;
} // namespace crosssection
} // namespace PROPOSAL

namespace PROPOSAL {
class CrossSectionDE2DX {
protected:
    size_t hash;
    std::shared_ptr<spdlog::logger> logger;
    double lower_energy_lim;

    CrossSectionDE2DX(size_t _hash);

public:
    CrossSectionDE2DX(crosssection::Parametrization<Medium> const&,
        ParticleDef const&, Medium const&, EnergyCutSettings const&,
        size_t hash = 0);

    CrossSectionDE2DX(crosssection::Parametrization<Component> const&,
        ParticleDef const&, Component const&, EnergyCutSettings const&,
        size_t hash = 0);

    virtual ~CrossSectionDE2DX() = default;

    virtual double Calculate(double energy) const = 0;

    virtual size_t GetHash() const noexcept { return hash; }
};

} // namespace PROPOSAL
