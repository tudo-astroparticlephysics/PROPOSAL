#pragma once
#include <memory>
#include "PROPOSAL/json_fwd.hpp"

namespace PROPOSAL {
    template <typename P, typename M>
    class CrossSection;
    struct ParticleDef;
    class Medium;
    class EnergyCutSettings;

    using cross_ptr = std::unique_ptr<CrossSection<ParticleDef, Medium>>;
}

namespace PROPOSAL {
    cross_ptr make_photonuclearQ2(const ParticleDef&, const Medium&,
                                  std::shared_ptr<const EnergyCutSettings>,
                                  bool, const nlohmann::json&);

    cross_ptr make_photonuclearreal(const ParticleDef&, const Medium&,
                                    std::shared_ptr<const EnergyCutSettings>,
                                    bool, const nlohmann::json&);
}