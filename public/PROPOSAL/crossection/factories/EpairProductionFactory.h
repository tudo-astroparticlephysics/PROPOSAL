
#pragma once


#include "PROPOSAL/methods.h"

namespace  PROPOSAL
{

class CrossSection;

class EpairProductionFactory
{
    public:

    struct Definition
    {
        Definition()
            : lpm_effect(true)
            , multiplier(1.0)
        {
        }

        bool lpm_effect;
        double multiplier;
    };

    // --------------------------------------------------------------------- //
    // Most general creation
    // --------------------------------------------------------------------- //

    CrossSection* CreateEpairProduction(const ParticleDef&,
                                        const Medium&,
                                        const EnergyCutSettings&,
                                        const Definition&) const;

    CrossSection* CreateEpairProduction(const ParticleDef&,
                                        const Medium&,
                                        const EnergyCutSettings&,
                                        const Definition&,
                                        InterpolationDef) const;

    // --------------------------------------------------------------------- //
    // Singleton pattern
    // --------------------------------------------------------------------- //

    static EpairProductionFactory& Get()
    {
        static EpairProductionFactory instance;
        return instance;
    }

    private:
    EpairProductionFactory();
    ~EpairProductionFactory();
};

} /*  PROPOSAL */

