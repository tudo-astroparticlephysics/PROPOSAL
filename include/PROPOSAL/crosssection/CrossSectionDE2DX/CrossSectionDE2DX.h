#pragma once

#include <cstddef>

namespace PROPOSAL {
namespace crosssection {
    template <typename Target> class Parametrization;
}
class ParticleDef;
class Medium;
class Component;
class EnergyCutSettings;
} // namespace PROPOSAL

namespace PROPOSAL {
class CrossSectionDE2DX {
protected:
    size_t hash;

public:
    CrossSectionDE2DX(crosssection::Parametrization<Component> const&,
        ParticleDef const&, Component const&, EnergyCutSettings const&);

    CrossSectionDE2DX(crosssection::Parametrization<Medium> const&,
        ParticleDef const&, Medium const&, EnergyCutSettings const&);

    virtual ~CrossSectionDE2DX() = default;

    virtual double Calculate(double energy) const = 0;

    virtual size_t GetHash() const noexcept { return hash; }
};

} // namespace PROPOSAL
