#pragma once

#include "PROPOSAL/secondaries/RegisteredInDefault.h"
#include "PROPOSAL/secondaries/annihilation/SingleDifferentialAnnihilation.h"

namespace PROPOSAL {
namespace secondaries {
    struct HeitlerAnnihilation
        : public secondaries::SingleDifferentialAnnihilation,
          public RegisteredInDefault<HeitlerAnnihilation> {

        HeitlerAnnihilation(const ParticleDef& p_def, const Medium& medium)
            : SingleDifferentialAnnihilation(crosssection::AnnihilationHeitler{}, p_def, medium, true)
        {
        }
    };
} // namespace secondaries
} // namespace PROPOSAL
