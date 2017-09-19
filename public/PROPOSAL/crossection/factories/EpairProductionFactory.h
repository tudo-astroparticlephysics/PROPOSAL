
#pragma once

#include <vector>
#include <map>
#include <string>

#include "PROPOSAL/crossection/parametrization/Parametrization.h"
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

    CrossSection* CreateEpairIntegral(const ParticleDef&,
                                      const Medium&,
                                      const EnergyCutSettings&,
                                      bool lpm,
                                      double multiplier) const;

    CrossSection* CreateEpairInterpolant(const ParticleDef&,
                                         const Medium&,
                                         const EnergyCutSettings&,
                                         bool lpm,
                                         double multiplier,
                                         InterpolationDef = InterpolationDef()) const;

    // --------------------------------------------------------------------- //
    // Most general creation
    // --------------------------------------------------------------------- //

    CrossSection* CreateEpairProduction(const ParticleDef&,
                                        const Medium&,
                                        const EnergyCutSettings&,
                                        bool lpm,
                                        double multiplier,
                                        bool interpolate,
                                        InterpolationDef = InterpolationDef()) const;

    CrossSection* CreateEpairProduction(const ParticleDef&,
                                        const Medium&,
                                        const EnergyCutSettings&,
                                        const Definition&,
                                        bool interpolate,
                                        InterpolationDef = InterpolationDef()) const;

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

