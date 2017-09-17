
#pragma once

#include <vector>
#include <map>
#include <string>

#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/geometry/GeometryFactory.h"

namespace  PROPOSAL
{

class UtilityFactory
{
    public:

    struct Definition : Utility::Definition
    {
        double e_cut;
        double v_cut;
        MediumFactory::Enum medium;
        double density_correction;

        Definition();
        ~Definition();
    };

    Utility* CreateUtility(const ParticleDef&, const Definition&);

    static UtilityFactory& Get()
    {
        static UtilityFactory instance;
        return instance;
    }

    private:
    UtilityFactory(){};
    ~UtilityFactory(){};
};

} /*  PROPOSAL */

