
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

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/particle/ParticleDef.h"

namespace PROPOSAL {

class CrossSection;
struct InterpolationDef;

typedef std::vector<std::shared_ptr<CrossSection>> Crosssections;

class Utility {
public:
    struct Definition {
        std::shared_ptr<InterpolationDef> inter_def;  // integration used
        std::shared_ptr<Scattering> scattering;       // no scattering
        std::shared_ptr<Crosssections> crosssections; // copy is expensive

        Definition(std::shared_ptr<Crosssections>,
            std::shared_ptr<Scattering> = nullptr,
            std::shared_ptr<InterpolationDef> = nullptr);
        ~Definition();
    };

    Utility(std::unique_ptr<Definition> utility_def)
        : utility_def(std::move(utility_def)){};

    DynamicData StochasticLoss(
        double particle_energy, std::array<double, 3> rnd);
    DynamicData AdvanceParticle(const DynamicData& p_condition,
        double final_energy, double sector_border);
    double DecayLength(const double initial_energy, const double rnd);
    double DecayEnergy(const double initial_energy, const double rnd);

private:
    Definition utility_def;

    std::unique_ptr<UtilityDecorator> displacement_calculator = nullptr;
    std::unique_ptr<UtilityDecorator> interaction_calculator = nullptr;
    std::unique_ptr<UtilityDecorator> decay_calculator = nullptr;
    std::unique_ptr<UtilityDecorator> cont_rand = nullptr;
    std::unique_ptr<UtilityDecorator> exact_time = nullptr;
};
} // namespace PROPOSAL

namespace PROPOSAL {
class UtilityDecorator {
    Crosssections crossections;

    public:
    UtilityDecorator(Crosssections cross) : crosssections(crosssections){};

    virtual double FunctionToIntegral(double energy) = 0;
    virtual double Calculate(double ei, double ef, double rnd) = 0;
    virtual double GetUpperLimit(double ei, double rnd) = 0;
};
} // namespace PROPOSAL
