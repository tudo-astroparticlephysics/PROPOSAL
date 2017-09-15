
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

class PropagationUtilityFactory
{
    public:

    struct Definition : PropagationUtility::Definition
    {
        double e_cut;
        double v_cut;
        bool do_interpolation;
        MediumFactory::Enum medium;
        double density_correction;

        Definition();
        ~Definition();
    };

    PropagationUtility* CreatePropagationUtility(const ParticleDef&, const Definition&);

    static PropagationUtilityFactory& Get()
    {
        static PropagationUtilityFactory instance;
        return instance;
    }

    private:
    PropagationUtilityFactory(){};
    ~PropagationUtilityFactory(){};
};

} /*  PROPOSAL */

