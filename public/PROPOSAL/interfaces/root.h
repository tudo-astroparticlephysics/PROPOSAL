#pragma once

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/Secondaries.h"

namespace PROPOSAL {
class Root {
    Root(std::string);

    void StoreDynamicData(const DynamicData& primary);
    void StoreSecondaries(const Secondaries& secondaries);
};
}
