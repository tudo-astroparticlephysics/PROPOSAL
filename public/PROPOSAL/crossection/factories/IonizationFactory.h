
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/


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

        bool operator==(const IonizationFactory::Definition& def) const
        {
            return multiplier == def.multiplier;
        }

        bool operator!=(const IonizationFactory::Definition& def) const
        {
            return !(*this == def);
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
