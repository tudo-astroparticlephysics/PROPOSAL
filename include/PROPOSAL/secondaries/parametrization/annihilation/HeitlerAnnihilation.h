#pragma once

#include "PROPOSAL/secondaries/parametrization/annihilation/SingleDifferentialAnnihilation.h"

namespace PROPOSAL {
namespace secondaries {
    struct HeitlerAnnihilation
        : public secondaries::SingleDifferentialAnnihilation,
          public DefaultSecondaries<HeitlerAnnihilation> {

        HeitlerAnnihilation(const ParticleDef& p_def, const Medium& medium)
            : SingleDifferentialAnnihilation(crosssection::AnnihilationHeitler{}, p_def, medium, true)
        {
        }
    };
} // namespace secondaries
} // namespace PROPOSAL
