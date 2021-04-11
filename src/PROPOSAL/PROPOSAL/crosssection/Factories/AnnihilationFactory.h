#pragma once
#include <nlohmann/json_fwd.hpp>
#include <memory>

namespace PROPOSAL {
struct CrossSectionBase;
struct ParticleDef;
class Medium;

using cross_ptr = std::unique_ptr<CrossSectionBase>;
}

namespace PROPOSAL {
cross_ptr make_annihilation(
    const ParticleDef&, const Medium&, bool, const nlohmann::json&);
}
