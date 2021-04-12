#pragma once
#include <nlohmann/json.hpp> // TODO: Use json_fwd.hpp as soon as https://github.com/conan-io/conan-center-index/pull/5149 is merged
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
