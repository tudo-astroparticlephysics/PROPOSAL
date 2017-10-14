
#pragma once


#include "PROPOSAL/methods.h"

namespace  PROPOSAL
{

class CrossSection;

class IonizationFactory
{
    public:

    struct Definition
    {
        Definition()
            : multiplier(1.0)
        {
        }

        double multiplier;
    };

    CrossSection* CreateIonizationIntegral(const ParticleDef&,
                                      const Medium&,
                                      const EnergyCutSettings&,
                                      double multiplier) const;

    CrossSection* CreateIonizationInterpolant(const ParticleDef&,
                                         const Medium&,
                                         const EnergyCutSettings&,
                                         double multiplier,
                                         InterpolationDef = InterpolationDef()) const;

    // --------------------------------------------------------------------- //
    // Most general creation
    // --------------------------------------------------------------------- //

    CrossSection* CreateIonization(const ParticleDef&,
                                        const Medium&,
                                        const EnergyCutSettings&,
                                        double multiplier,
                                        bool interpolate,
                                        InterpolationDef = InterpolationDef()) const;

    CrossSection* CreateIonization(const ParticleDef&,
                                        const Medium&,
                                        const EnergyCutSettings&,
                                        const Definition&,
                                        bool interpolate,
                                        InterpolationDef = InterpolationDef()) const;

    // --------------------------------------------------------------------- //
    // Singleton pattern
    // --------------------------------------------------------------------- //

    static IonizationFactory& Get()
    {
        static IonizationFactory instance;
        return instance;
    }

    private:
    IonizationFactory();
    ~IonizationFactory();
};

} /*  PROPOSAL */

