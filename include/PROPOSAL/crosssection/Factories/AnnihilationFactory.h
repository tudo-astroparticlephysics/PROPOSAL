#pragma once
#include <memory>
#include "PROPOSAL/json_fwd.hpp"

namespace PROPOSAL {
    template <typename P, typename M>
    class CrossSection;
    struct ParticleDef;
    class Medium;

    using cross_ptr = std::unique_ptr<CrossSection<ParticleDef, Medium>>;
}

namespace PROPOSAL {
    cross_ptr make_annihilation(const ParticleDef&, const Medium&, bool,
                                const nlohmann::json&);
}