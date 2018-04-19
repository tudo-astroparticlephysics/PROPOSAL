
#pragma once

#include "PROPOSAL/methods.h"

namespace PROPOSAL {

class CrossSection;
struct ParticleDef;
class Medium;
class EnergyCutSettings;

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

    // --------------------------------------------------------------------- //
    // Most general creation
    // --------------------------------------------------------------------- //

    CrossSection* CreateIonization(const ParticleDef&,
                                   const Medium&,
                                   const EnergyCutSettings&,
                                   const Definition&) const;

    CrossSection* CreateIonization(const ParticleDef&,
                                   const Medium&,
                                   const EnergyCutSettings&,
                                   const Definition&,
                                   InterpolationDef) const;

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

} // namespace PROPOSAL
