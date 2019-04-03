
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

#include <vector>

#include "PROPOSAL/crossection/factories/BremsstrahlungFactory.h"
#include "PROPOSAL/crossection/factories/EpairProductionFactory.h"
#include "PROPOSAL/crossection/factories/IonizationFactory.h"
#include "PROPOSAL/crossection/factories/PhotonuclearFactory.h"
#include "PROPOSAL/crossection/factories/MupairProductionFactory.h"
#include "PROPOSAL/crossection/factories/WeakInteractionFactory.h"

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/particle/ParticleDef.h"

namespace PROPOSAL {

class CrossSection;
class Medium;

struct InterpolationDef;

class Utility
{
public:
    struct Definition
    {
        BremsstrahlungFactory::Definition brems_def;
        PhotonuclearFactory::Definition photo_def;
        EpairProductionFactory::Definition epair_def;
        IonizationFactory::Definition ioniz_def;
        MupairProductionFactory::Definition mupair_def;
        WeakInteractionFactory::Definition weak_def;


        bool operator==(const Utility::Definition& utility_def) const;
        bool operator!=(const Utility::Definition& utility_def) const;

        Definition();
        ~Definition();
    };

public:
    Utility(const ParticleDef&, const Medium&, const EnergyCutSettings&, Definition);
    Utility(const ParticleDef&, const Medium&, const EnergyCutSettings&, Definition, const InterpolationDef&);
    Utility(const std::vector<CrossSection*>&);
    Utility(const Utility&);
    virtual ~Utility();

    bool operator==(const Utility& utility) const;
    bool operator!=(const Utility& utility) const;

    const ParticleDef& GetParticleDef() const { return particle_def_; }
    const Medium& GetMedium() const { return *medium_; }
    const std::vector<CrossSection*>& GetCrosssections() const { return crosssections_; }

protected:
    Utility& operator=(const Utility&); // Undefined & not allowed

    // --------------------------------------------------------------------- //
    // Protected members
    // --------------------------------------------------------------------- //

    ParticleDef particle_def_;
    Medium* medium_;
    EnergyCutSettings cut_settings_;

    std::vector<CrossSection*> crosssections_;
};

class UtilityDecorator
{
public:
    UtilityDecorator(const Utility&);

    // Copy constructors
    UtilityDecorator(const UtilityDecorator&);
    virtual UtilityDecorator* clone(const Utility&) const = 0;

    virtual ~UtilityDecorator();

    bool operator==(const UtilityDecorator& utility_decorator) const;
    bool operator!=(const UtilityDecorator& utility_decorator) const;

    // Methods
    virtual double FunctionToIntegral(double energy);
    virtual double Calculate(double ei, double ef, double rnd) = 0;
    virtual double GetUpperLimit(double ei, double rnd)        = 0;

    const Utility& GetUtility() const { return utility_; }

protected:
    UtilityDecorator& operator=(const UtilityDecorator&); // Undefined & not allowed

    // Implemented in child classes to be able to use equality operator
    virtual bool compare(const UtilityDecorator&) const = 0;

    const Utility& utility_;
};

} // namespace PROPOSAL
