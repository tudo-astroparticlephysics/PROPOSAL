#pragma once
#include <memory>
#include <nlohmann/json_fwd.hpp>

namespace PROPOSAL {
    struct CrossSectionBase;
    struct ParticleDef;
    class Medium;
    class EnergyCutSettings;

    using cross_ptr = std::unique_ptr<CrossSectionBase>;
}

namespace PROPOSAL {
    cross_ptr make_photonuclearQ2(const ParticleDef&, const Medium&,
                                  std::shared_ptr<const EnergyCutSettings>,
                                  bool, const nlohmann::json&);

    cross_ptr make_photonuclearreal(const ParticleDef&, const Medium&,
                                    std::shared_ptr<const EnergyCutSettings>,
                                    bool, const nlohmann::json&);
}
