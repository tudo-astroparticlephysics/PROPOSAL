
#pragma once


#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/medium/MediumFactory.h"

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

