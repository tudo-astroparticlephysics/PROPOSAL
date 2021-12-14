#pragma once
#include <memory>
#include <nlohmann/json_fwd.hpp>

namespace PROPOSAL {
    struct CrossSectionBase;
    struct ParticleDef;
    class Medium;

    using cross_ptr = std::unique_ptr<CrossSectionBase>;
}

namespace PROPOSAL {
    cross_ptr make_photoproduction(const ParticleDef&, const Medium&,
                                   const nlohmann::json&);
}
