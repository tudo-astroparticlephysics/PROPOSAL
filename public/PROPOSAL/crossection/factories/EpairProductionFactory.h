
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
class ParticleDef;
class Medium;
class EnergyCutSettings;

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

} // namespace PROPOSAL
